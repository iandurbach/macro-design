library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(mapview)
library(leaflet)
library(htmltools)
library(raster)
library(secr)
library(secrdesign)
library(purrr)

source("cv-full-sim/simulate_capthists.R")
source("cv-full-sim/cv-utils.R")

# mean function for non-euc distance calcs
mymean = function(x) exp(mean(log(x))) 

my_crs <- "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs"

#########################
### reading in single mesh and trap array
#########################

path_to_mesh <- "data/TostMask.Rds"
path_to_traps <- "data/TostTraps.Rds"

nsims <- 2
mesh_reduce <- 3
n_pts <- 20
dens_per_100km2 <- NULL
b_ac <- 0
cov_ac <- c("stdGC")
b_con <- 0
cov_con <- c("stdGC")
lambda0 <- c(2)
sigma <- c(4000)
n_occasions <- c(1)

mesh <- readRDS(path_to_mesh)
traps <- readRDS(path_to_traps)

#########################
### reading in multiple meshes and trap arrays
#########################

path_to_mesh_3regions <- "data/allmasks.Rds"
path_to_traps_3regions <- "data/alltraps.Rds"

mesh <- readRDS(path_to_mesh_3regions)
traps <- readRDS(path_to_traps_3regions)

nsims <- 3
mesh_reduce <- c(3,3,3)
n_pts <- c(NULL,NULL,NULL)
dens_per_100km2 <- c(1.5, 1.5, 1.5)
b_ac <- c(0,0,0)
cov_ac <- c("stdGC", "stdGC", "stdGC")
b_con <- c(0,0,0) # positive means higher values more conductive
cov_con <- c("stdGC", "stdGC", "stdGC")
lambda0 <- c(1.5,1.5,1.5)
sigma <- c(8000,8000,8000)
n_occasions <- c(1,1,1)

#########################
## end of reading in data
#########################
#########################
#########################

if(!is.null(names(mesh))){
  mesh <- list(mesh)
  traps <- list(traps)
}

# reduce mesh size (set mesh_reduce = rep(1, nregions) for no reduction)
mesh <- reduce_mesh_size(mesh = mesh, mesh_reduce = mesh_reduce, cov_ac = cov_ac, cov_con = cov_con)

mesh_sf <- list()
traps_sf <- list()
for(i in 1:length(mesh)){
  # make the mesh an sf object
  mesh_tmp <- raster::raster(mesh[[i]], covariate = cov_ac[i], values = 1, crs = NA)
  mesh_tmp <- as(mesh_tmp,'SpatialPolygonsDataFrame')
  mesh_sf[[i]] <- st_as_sf(mesh_tmp, crs = my_crs)
  mesh_sf[[i]] <- mesh_sf[[i]] %>% st_set_crs(my_crs)
  rm(mesh_tmp)
  
  # make the traps an sf object
  traps_sf[[i]] <- st_as_sf(traps[[i]], coords = c("x","y"))
  traps_sf[[i]] <- traps_sf[[i]] %>% st_set_crs(my_crs)
  
}

m <- leaflet() %>% 
  addTiles() %>%
  #  flyToBounds(85, 30, 86, 53) %>%
  addPolygons(data = st_transform(purrr::reduce(mesh_sf, rbind), crs = 4326), fillOpacity = 0.2, weight = 1) %>%
  addCircleMarkers(data = st_transform(purrr::reduce(traps_sf, rbind), crs = 4326), radius = 1, color = "red") 

m

# simulate several capture histories in each region

# turn the matrix of trap co-ordinates into a trap object
for(i in 1:length(mesh)){
  traps[[i]] <- read.traps(data = traps[[i]], detector = "count", binary.usage = FALSE)
}

sim_ch <- purrr::map(1:nsims, simulate_capthists, mesh = mesh, traps = traps, n_pts = n_pts, dens_per_100km2 = dens_per_100km2,
                     b_ac = b_ac, cov_ac = cov_ac, b_con = b_con,
                     cov_con = cov_con, lambda0 = lambda0, sigma = sigma, n_occasions = n_occasions,
                     myLCdist = mymean)

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
summary_sim_ch

# secr.fit to each region separately
mods <- purrr::map(sim_ch, ~ pmap(list(capthist = .$ch, mask = mesh, start = .$startvals),
                                                      secr.fit, detectfn = "HHN", model = list(D ~ 1)))

# furrr not stable with RStudio
# future::plan(multiprocess)
# mods <- future_map(sim_ch, ~ future_pmap(list(capthist = .$ch, mask = mesh, start = .$startvals), 
#                                          secr.fit, detectfn = "HHN", model = list(D ~ 1)), .progress = TRUE)

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

# make cv into nice table
cv_res <- data.frame(sim_id = rep(1:nrow(expected_cv), times = ncol(expected_cv)),
                     region = rep(1:ncol(expected_cv), each = nrow(expected_cv)),
                     E_cv = as.vector(expected_cv),
                     R_cv = as.vector(realised_cv))
cv_table <- cv_res %>% group_by(region) %>% summarize(CV_E.N = mean(E_cv), CV_R.N = mean(R_cv))
globcv_table <- data.frame(region = "all", CV_E.N = mean(global_expected_cv), CV_R.N = mean(global_realised_cv))

rbind(cv_table, globcv_table)

