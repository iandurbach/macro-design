# This app calculates approximate CV values for different numbers of sites and cameras, using
# the approximations in Efford and Boulanger (2019).
# The approximation assumes constant animal density and likely underestimates CV when 
# density is non-uniform. The app includes a 'correction factor' that penalizes the CV for extra-Poisson 
# variation in density but this is an experimental feature and it is not clear what values are most appropriate. 
# This app is intended for getting a rough idea of appropriate numbers of sites and cameras, and not for the 
# selection of a final design. This should always be assessed by more comprehensive simulations. A companion
# app is available for this purpose.

library(dplyr)
library(ggplot2)
library(plotly)
library(secrdesign)
library(viridis)

range_nsurveys <- c(1,10)
range_ntraps <- c(20,40)

sigma <- 4000
lambda0 <- 0.5
D_per_100km2 <- 0.25

D <- D_per_100km2 / 10000 # turn D per 100km2 into D per ha

trapspacing_mult <- 1 # traps spaced a multiple of sigma apart
CF <- 1

poss_n_hal <- round(seq(from = range_nsurveys[1], to = range_nsurveys[2], length.out = min(8, diff(range_nsurveys) + 1)), 0)
poss_n_cam <- round(seq(from = range_ntraps[1], to = range_ntraps[2], length.out = min(8, diff(range_ntraps) + 1)), 0)

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
res$CV <- res$rotRSE * CF

limmin <- min(0, min(res$CV))
limmax <- max(1, max(res$CV))

p <- res %>%
  mutate(txt = paste("CV:", CV, "\n", "E(n):", En, "\n", "E(r):", Er)) %>%
  ggplot(aes(x = nrepeats, y = ntraps)) +
  geom_raster(aes(fill = CV, text = txt)) +
  xlab("Number of sites") +
  ylab("Number of cameras per site") +
  scale_fill_viridis(limits = c(limmin,limmax), direction = -1) +
  theme(legend.position = "bottom")

if((min(res$CV) < 0.3) & (max(res$CV) > 0.3)){
  p <- p + geom_contour(aes(z = CV), breaks = c(0.3), colour = "white")
}

if((min(res$CV) < 0.2) & (max(res$CV) > 0.2)){
  p <- p + geom_contour(aes(z = CV), breaks = c(0.2), colour = "gray80")
}

ggplotly(p, tooltip = "txt")

