## Macro-design app -- generates Halton points over a defined study area, with
## stratification

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(sf)
library(mapview)
library(leaflet)
library(BalancedSampling)

source("generate-macro-sites/generate_Halton_pts.R")
source("generate-macro-sites/macro-utils.R")
source("generate-macro-sites/macro-design-fn.R")

macro_design(path_to_study_area = "data/Occu15x15_PotSL_wstrata.shp", path_to_surveyed_pts = NULL, buffer = 10000,
             strat_var = "RASTERVALU", strat_breaks = c(0,0.5,0.8,1), strat_proportions = NULL,
             initial_seed = 1234, n_existing_pts = 0, n_new_pts = 50,
             basemap = "Esri.WorldStreetMap")
