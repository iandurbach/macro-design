simulate_capthists <- function(mesh, traps, n_pts, dens_per_100km2, b_ac, cov_ac, b_con, cov_con, lambda0, sigma, n_occasions, ...){
  
  ## generate non-uniform AC density
  
  # ac density ~ ruggedness, so add this as a mesh covariate.
  # b_ac is a parameter controlling strength of relationship
  simulated_points_Dcov <- list()
  Dcov_for_sim <- list()
  ch <- list()
  n <- list()
  r <- list()
  startvals <- list()
  for(i in 1:length(mesh)){
    covariates(mesh[[i]])$Dac <- exp(b_ac[i] * covariates(mesh[[i]])[[cov_ac[i]]])
    if(!is.null(n_pts[i])){
      Dcov_for_sim[[i]] <- n_pts[i] / attr(mesh[[i]], "area") * (covariates(mesh[[i]])$Dac / sum(covariates(mesh[[i]])$Dac)) 
    }else{
      implied_n_pts <- (dens_per_100km2[i] / 10000) * attr(mesh[[i]], "area") * length(mesh[[i]]$x)
      Dcov_for_sim[[i]] <- implied_n_pts / attr(mesh[[i]], "area") * (covariates(mesh[[i]])$Dac / sum(covariates(mesh[[i]])$Dac))
    }
    simulated_points_Dcov[[i]] <- sim.popn(D = Dcov_for_sim[[i]], 
                                           core = mesh[[i]], 
                                           model2D = "IHP",
                                           Ndist = "fixed")
    
    ## generate non-Euclidean conductance
    
    # conductance ~ ruggedness, so add this as a mesh covariate.
    # create pixel-specific cost/friction and assign to the simulated popn objects 
    # b_con is a parameter controlling strength of relationship
    # make sure you have specified "my_mean" function correctly!
    
    covariates(mesh[[i]])$noneuc <- exp(b_con[i] * covariates(mesh[[i]])[[cov_con[i]]])
    attr(simulated_points_Dcov[[i]], "mask") <- mesh[[i]]
    
    ## simulate a capture history 
    ch[[i]] <- sim.capthist(traps = traps[[i]], pop = simulated_points_Dcov[[i]], 
                            userdist = myLCdist, 
                            noccasions = n_occasions[i],
                            detectpar = list(lambda0 = lambda0[i], sigma = sigma[i]), 
                            detectfn = "HHN")
    
    # make starting values for secr.fit
    startvals[[i]] <- list(D = mean(Dcov_for_sim[[i]]), lambda0 = lambda0[i], sigma = sigma[i])
    
  }
  
  return(list(ch = ch, startvals = startvals))
  
}