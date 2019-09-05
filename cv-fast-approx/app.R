library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(secrdesign)
library(viridis)

ui <- fluidPage(

  titlePanel("Fast CV approximations for designing camera trap surveys"),

  sidebarLayout(
    sidebarPanel(
      sliderInput("nsurveys",
                  "Number of sites",
                  min = 1,
                  max = 50,
                  value = c(1,10)),
      sliderInput(inputId = "ntraps",
                   label = "Number of cameras per site",
                  min = 1,
                  max = 80,
                   value = c(20,40)),
      numericInput(inputId = "D",
                  label = "Expected D / 100km2",
                  min = 0.01,
                  step = 0.05,
                  value = c(0.25)),
      numericInput(inputId = "lambda0",
                   label = "Expected lambda",
                   min = 0.001,
                   step = 0.001,
                   value = 0.5),
      numericInput(inputId = "sigma",
                   label = "Expected sigma",
                   min = 0.001,
                   step = 100,
                   value = 4000),
      numericInput(inputId = "trapspacing",
                   label = "Trap spacing (multiple of sigma)",
                   min = 0.001,
                   step = 0.1,
                   value = 1),
      sliderInput(inputId = "CF",
                   label = "Correction factor (experimental)",
                  min = 1,
                  max = 2,
                  step = 0.1,
                  value = 1),
      actionButton(inputId = "calc", label = "Calculate CV")
    ),

    mainPanel(
      p("This app calculates approximate CV values for different numbers of sites and cameras, using
      the approximations in Efford and Boulanger (2019)."),
      p("The approximation assumes constant animal density and likely underestimates CV when 
      density is non-uniform. The app includes a 'correction factor' that penalizes the CV for extra-Poisson 
        variation in density but this is an experimental feature and it is not clear what values are most appropriate."),
      p("This app is intended for getting a rough idea of appropriate numbers of sites and cameras, and not for the 
      selection of a final design. This should always be assessed by more comprehensive simulations. A companion
        app is available for this purpose."),
      plotlyOutput("cv_res")
    )
  )
)

server <- function(input, output) {

  cv_calcs <- eventReactive(input$calc, {

    sigma <- input$sigma
    lambda0 <- input$lambda0
    D <- input$D / 10000 # turn D per 100km2 into D per ha

    trapspacing_mult <- input$trapspacing # traps spaced a multiple of sigma apart

    poss_n_hal <- round(seq(from = input$nsurveys[1], to = input$nsurveys[2], length.out = min(8, diff(input$nsurveys) + 1)), 0)
    poss_n_cam <- round(seq(from = input$ntraps[1], to = input$ntraps[2], length.out = min(8, diff(input$ntraps) + 1)), 0)

    mydetectors <- list()

    for(i in 1:length(poss_n_cam)){
      n_cam <- poss_n_cam[i]
      n_cam_x <- ceiling(sqrt(n_cam))
      mydetectors[[i]] <- read.traps(data = data.frame(expand.grid(x = seq(from = 0, by = sigma * trapspacing_mult, length.out = n_cam_x),
                                                                   y = seq(from = 0, by = sigma * trapspacing_mult, length.out = n_cam_x))[1:n_cam,],
                                                       row.names = 1:n_cam), detector = "count")
    }

    scen <- make.scenarios(trapsindex = 1:length(poss_n_cam), D = D, sigma = sigma, nrepeats = poss_n_hal, detectfn = "HHN", noccasions = 1, lambda0 = lambda0)

    res <- scenarioSummary(scen, mydetectors)
    res$ntraps <- poss_n_cam[res$trapsindex]
    res$CV <- res$rotRSE * input$CF

    return(list(res = res, scen = scen))

  })

  output$cv_res <- renderPlotly({

    limmin <- min(0, min(cv_calcs()$res$CV))
    limmax <- max(1, max(cv_calcs()$res$CV))

    p <- cv_calcs()$res %>%
      mutate(txt = paste("CV:", CV, "\n", "E(n):", En, "\n", "E(r):", Er)) %>%
      ggplot(aes(x = nrepeats, y = ntraps)) +
      geom_raster(aes(fill = CV, text = txt)) +
      xlab("Number of sites") +
      ylab("Number of cameras per site") +
      scale_fill_viridis(limits = c(limmin,limmax), direction = -1) +
      theme(legend.position = "bottom")

    if((min(cv_calcs()$res$CV) < 0.3) & (max(cv_calcs()$res$CV) > 0.3)){
      p <- p + geom_contour(aes(z = CV), breaks = c(0.3), colour = "red")
    }

    if((min(cv_calcs()$res$CV) < 0.2) & (max(cv_calcs()$res$CV) > 0.2)){
      p <- p + geom_contour(aes(z = CV), breaks = c(0.2), colour = "blue")
    }

    ggplotly(p, tooltip = "txt")

  })

}

shinyApp(ui = ui, server = server)
