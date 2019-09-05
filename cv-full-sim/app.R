library(shiny)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(sf)
library(mapview)
library(leaflet)
library(htmltools)
library(raster)
library(secr)
library(secrdesign)

source("simulate_capthists.R")
source("cv-utils.R")

ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("map", width = "100%", height = "100%"),
  absolutePanel(top = 100, left = 10,
                fluidRow(
                  column(3, fileInput("trap_file",
                                      "File for traps")),
                  column(3, fileInput("region_file",
                                      "File for meshes"))),
                fluidRow(
                  column(3, numericInput("nsims",
                                         "Simulations to run",
                                         min = 1,
                                         value = 3)),
                  column(3, numericInput(inputId = "dens_per_100km2",
                                         label = "Expected D/100km2",
                                         min = 0,
                                         value = 1.5))),
                fluidRow(
                  column(3, textInput(inputId = "cov_ac",
                                      label = "Name of density covariate",
                                      value = "stdGC")),
                  column(3, numericInput(inputId = "b_ac",
                                         label = "Coeff. of density covariate",
                                         min = 0,
                                         value = 0))),
                fluidRow(
                  column(3, numericInput(inputId = "lambda0",
                                         label = "Expected lambda",
                                         min = 0,
                                         value = 1.5)),
                  column(3, numericInput(inputId = "sigma",
                                         label = "Expected sigma",
                                         min = 0,
                                         value = 4000))),
                fluidRow(
                  column(3, textInput(inputId = "cov_con",
                                      label = "Name of conductance covariate",
                                      value = "stdGC")),
                  column(3, numericInput(inputId = "b_con",
                                         label = "Coeff. of conductance covariate",
                                         min = 0,
                                         value = 0))),
                selectInput("basemap","Base map",
                            choices = list("Esri.WorldStreetMap", "Esri.WorldImagery", "OpenStreetMap"),
                            selected = "Esri.WorldStreetMap"),
                div(actionButton(inputId = "add_region", label = "Add region"),
                    actionButton(inputId = "run_sim_ch", label = "Simulate histories"),
                    actionButton(inputId = "sim_mods", label = "Run SCR"),
                    style = "margin:2px"),
  checkboxInput(inputId = "reduce_mesh", label = "Reduce mesh size")),
  absolutePanel(top = 10, right = 10, 
                fluidRow(
                  column(6, tableOutput("ch")))),
  absolutePanel(top = 200, right = 10, 
                fluidRow(
                  column(12, tableOutput("cv")))),
  absolutePanel(bottom = 20, left = 10,
                div(downloadButton("downloadHists", "Download capthists"),
                    downloadButton("downloadModels", "Download all")))
)

# Define server logic
server <- function(input, output, session) {
  
  all_region_mesh <- list()
  all_region_traps <- list()
  all_region_pars <- data.frame(mesh_reduce = as.numeric(),
                                dens_per_100km2 = as.numeric(),
                                b_ac = as.numeric(),
                                cov_ac = as.character(),
                                b_con = as.numeric(),
                                cov_con = as.character(),
                                lambda0 = as.numeric(),
                                sigma = as.numeric(),
                                n_occasions = as.numeric())
  
  region_data <- eventReactive(input$add_region, {
    
    # work out what UTM zone we're in based on the centre of the visible map
    my_utm_zone <- floor((input$map_center$lng + 180) / 6) + 1
    my_crs <- paste0("+proj=utm +zone=", my_utm_zone, " +datum=WGS84 +units=m +no_defs")
    
    region_mesh <- readRDS(input$region_file$datapath)
    region_traps <- readRDS(input$trap_file$datapath)
    
    if(!is.null(names(region_mesh))){
      region_mesh <- list(region_mesh)
      region_traps <- list(region_traps)
    }
    
    # reduce mesh size (set mesh_reduce = rep(1, nregions) for no reduction)
    if(input$reduce_mesh){
      mesh_reduce_factor <- rep(3, length(region_mesh))
    } else { mesh_reduce_factor <- rep(1, length(region_mesh)) }

    region_mesh <- reduce_mesh_size(mesh = region_mesh, mesh_reduce = mesh_reduce_factor, 
                                    cov_ac = rep(input$cov_ac, length(region_mesh)), 
                                    cov_con = rep(input$cov_con, length(region_mesh)))
    
    region_mesh_sf <- list()
    region_traps_sf <- list()
    for(i in 1:length(region_mesh)){
      # make the mesh an sf object
      region_mesh_tmp <- raster::raster(region_mesh[[i]], covariate = input$cov_ac, values = 1, crs = NA)
      region_mesh_tmp <- as(region_mesh_tmp,'SpatialPolygonsDataFrame')
      region_mesh_sf[[i]] <- st_as_sf(region_mesh_tmp, crs = my_crs)
      region_mesh_sf[[i]] <- region_mesh_sf[[i]] %>% st_set_crs(my_crs)
      rm(region_mesh_tmp)
      
      # make the traps an sf object
      region_traps_sf[[i]] <- st_as_sf(region_traps[[i]], coords = c("x","y"))
      region_traps_sf[[i]] <- region_traps_sf[[i]] %>% st_set_crs(my_crs)
    }
    
    region_pars <- data.frame(mesh_reduce = rep(1,length(region_mesh)),
                              dens_per_100km2 = rep(input$dens_per_100km2,length(region_mesh)),
                              b_ac = rep(input$b_ac,length(region_mesh)),
                              cov_ac = rep(input$cov_ac,length(region_mesh)),
                              b_con = rep(input$b_con,length(region_mesh)),
                              cov_con = rep(input$cov_con,length(region_mesh)),
                              lambda0 = rep(input$lambda0,length(region_mesh)),
                              sigma = rep(input$sigma,length(region_mesh)),
                              n_occasions = rep(1,length(region_mesh)))
    
    return(list(region_mesh = region_mesh, region_mesh_sf = region_mesh_sf, 
                region_traps = region_traps, region_traps_sf = region_traps_sf,
                region_pars = region_pars))
    
  })
  
  update_region_lists <- reactiveValues()
  
  observeEvent(region_data(), {
    update_region_lists$all_region_mesh <- c(isolate(update_region_lists$all_region_mesh), isolate(region_data()$region_mesh))
    update_region_lists$all_region_traps <- c(isolate(update_region_lists$all_region_traps), isolate(region_data()$region_traps))
    update_region_lists$all_region_pars <- rbind(isolate(update_region_lists$all_region_pars), isolate(region_data()$region_pars))
  })
  
  capture_histories <- eventReactive(input$run_sim_ch, {
    
    traps <- list()
    for(i in 1:length(update_region_lists$all_region_mesh)){
      traps[[i]] <- read.traps(data = update_region_lists$all_region_traps[[i]], detector = "count", binary.usage = FALSE)
    }
    
    # simulate several capture histories in each region
    sim_ch <- purrr::map(1:input$nsims, simulate_capthists,
                         mesh = update_region_lists$all_region_mesh,
                         traps = traps,
                         n_pts = rep(NULL, length(traps)),
                         dens_per_100km2 = update_region_lists$all_region_pars$dens_per_100km2,
                         b_ac = update_region_lists$all_region_pars$b_ac,
                         cov_ac = update_region_lists$all_region_pars$cov_ac,
                         b_con = update_region_lists$all_region_pars$b_con,
                         cov_con = update_region_lists$all_region_pars$cov_con,
                         lambda0 = update_region_lists$all_region_pars$lambda0,
                         sigma = update_region_lists$all_region_pars$sigma,
                         n_occasions = rep(1, length(traps)))
    
    
    # inspect simulated capthists
    chs <- purrr::map(sim_ch, ~.x[["ch"]])
    n_detections <- c()
    n_recaps <- c()
    for(i in 1:length(chs)){
      for(j in 1:length(chs[[1]])){
        n_dets <- c(n_detections, summary(chs[[i]][[j]])$counts$Total[4])
        n_recaps <- c(n_recaps, summary(chs[[i]][[j]])$counts$Total[6])
      }
    }
    summary_sim_ch <- data.frame(sim_id = rep(1:length(chs), each = length(chs[[1]])),
                                 region = rep(1:length(chs[[1]]), times = length(chs)),
                                 n_dets, n_recaps)
    
    return(list(sim_ch = sim_ch, summary_sim_ch = summary_sim_ch))
    
  })
  
  fitted_models <- eventReactive(input$sim_mods, {
    
    mods <- purrr::map(capture_histories()$sim_ch, ~ pmap(list(capthist = .$ch, mask = update_region_lists$all_region_mesh, start = .$startvals),
                                                          secr.fit, detectfn = "HHN", model = list(D ~ 1)), .progress = TRUE)
    
    # calculate cv's
    
    # calculate Nhat and se.Nhat for each fitted model
    expected_NpR_est <- matrix(NA, nrow = length(mods), ncol = length(mods[[1]]))
    expected_NpR_se <- matrix(NA, nrow = length(mods), ncol = length(mods[[1]]))
    realised_NpR_est <- matrix(NA, nrow = length(mods), ncol = length(mods[[1]]))
    realised_NpR_se <- matrix(NA, nrow = length(mods), ncol = length(mods[[1]]))
    for(i in 1:length(mods)){
      N_per_region <- purrr::map(mods[[i]], region.N)
      expected_NpR_est[i, ] <- N_per_region %>% purrr::map(~.x[["estimate"]][1]) %>% unlist()
      expected_NpR_se[i, ] <- N_per_region %>% purrr::map(~.x[["SE.estimate"]][1]) %>% unlist()
      realised_NpR_est[i, ] <- N_per_region %>% purrr::map(~.x[["estimate"]][2]) %>% unlist()
      realised_NpR_se[i, ] <- N_per_region %>% purrr::map(~.x[["SE.estimate"]][2]) %>% unlist()
    }
    
    # calculate CV of each Nhat i.e. each region, each simulation
    expected_cv <- expected_NpR_se / expected_NpR_est
    realised_cv <- realised_NpR_se / realised_NpR_est
    
    # calculate overall expected Nhat, Nhat.se, and cv for each simulation
    global_expected_N <- apply(expected_NpR_est, 1, sum)
    global_expected_N_se <- apply(expected_NpR_se, 1, function(x) sqrt(sum(x^2)))
    global_expected_cv <- global_expected_N_se / global_expected_N
    
    # calculate overall realised Nhat, Nhat.se, and cv for each simulation
    global_realised_N <- apply(realised_NpR_est, 1, sum)
    global_realised_N_se <- apply(realised_NpR_se, 1, function(x) sqrt(sum(x^2)))
    global_realised_cv <- global_realised_N_se / global_realised_N
    
    globalcv <- data.frame(CV_E.N = global_expected_cv, CV_R.N = global_realised_cv)
    
    # make cv into df
    cv_res <- data.frame(sim_id = rep(1:nrow(expected_cv), times = ncol(expected_cv)),
                         region = rep(1:ncol(expected_cv), each = nrow(expected_cv)),
                         E_cv = as.vector(expected_cv),
                         R_cv = as.vector(realised_cv))
    
    return(list(mods = mods, cv = cv_res, globalcv = globalcv))
    
  })
  
  # only include aspects of the map that won't need to change dynamically
  output$map <- renderLeaflet({
    leaflet() %>% addTiles() %>%
      flyToBounds(100, 42, 103, 50)
  })
  
  # purrr::reduce(mesh_sf, rbind)
  
  observe({
    leafletProxy("map") %>%
      clearTiles() %>%
      addTiles() %>%
      addProviderTiles(eval(parse(text=paste0("providers$",input$basemap))))
  })
  
  observe({
    leafletProxy("map") %>%
      addPolygons(data = st_transform(purrr::reduce(region_data()$region_mesh_sf, rbind), crs = 4326), fillOpacity = 0.2, weight = 1) %>%
      addCircleMarkers(data = st_transform(purrr::reduce(region_data()$region_traps_sf, rbind), crs = 4326), radius = 1, color = "red")
  })
  
  output$ch <- renderTable({
    capture_histories()$summary_sim_ch %>% group_by(region) %>% summarize(N.min = min(n_dets),
                                                                          N.mean = mean(n_dets),
                                                                          N.max = max(n_dets),
                                                                          R.min = min(n_recaps),
                                                                          R.mean = mean(n_recaps),
                                                                          R.max = max(n_recaps))
  }, digits = 0)
  
  output$cv <- renderTable({
    rbind(fitted_models()$cv %>% group_by(region) %>% summarize(CV_E.N = mean(E_cv), CV_R.N = mean(R_cv)), 
          data.frame(region = "all", CV_E.N = mean(fitted_models()$globalcv$CV_E.N),
                     CV_R.N = mean(fitted_models()$globalcv$CV_R.N)))
  })
  
  # downloadable RData file with capture histories
  output$downloadHists <- downloadHandler(
    filename = function() {
      paste0("simulated_capthists.RData")
    },
    content = function(file) {
      
      mesh <- update_region_lists$all_region_mesh
      traps <- update_region_lists$all_region_traps
      pars <- update_region_lists$all_region_pars
      capthists <- capture_histories()$sim_ch

      save(mesh, traps, pars, capthists, file = file)
      }
  )
  
  output$downloadModels <- downloadHandler(
    filename = function() {
      paste0("simulated_models.RData")
    },
    content = function(file) {
      
      mesh <- update_region_lists$all_region_mesh
      traps <- update_region_lists$all_region_traps
      pars <- update_region_lists$all_region_pars
      capthists <- capture_histories()$sim_ch
      models <- fitted_models()$mods
      cvs <- fitted_models()$cv
      globalcv <- fitted_models()$globalcv
      
      save(mesh, traps, pars, capthists, models, cvs, globalcv, file = file)
    }
  )
  
  
}

# Run the application
shinyApp(ui = ui, server = server)

