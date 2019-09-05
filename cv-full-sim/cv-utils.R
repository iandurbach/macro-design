### Utility function for simulating-cv.

# mean function for non-euc distance calcs
mymean <- function(x) exp(mean(log(x))) 

# least cost distance function
myLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = mymean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

# function for reducing the number of mesh points by a factor mesh_reduce
reduce_mesh_size <- function(mesh, mesh_reduce, cov_ac, cov_con){
  
  # reduce mesh by some factor for faster processing (leave at 1 for as is)
  red_mesh <- list()
  for(i in 1:length(mesh)){
    redmesh_cov_ac <- secr::raster(mesh[[i]], cov_ac[i])
    redmesh_cov_ac <- raster::aggregate(redmesh_cov_ac, fact = mesh_reduce[i], fun = mean)
    ac_name = quo_name(cov_ac[i])
    redmesh_cov_ac_df <- as.data.frame(redmesh_cov_ac) %>% dplyr::rename(!!ac_name := layer)
    red_df <- cbind(as.data.frame(coordinates(redmesh_cov_ac)), redmesh_cov_ac_df, D = 1) 
    red_df <- red_df %>% dplyr::filter_at(cov_ac[i], all_vars(!is.na(.)))
    
    if(cov_con[i] != cov_ac[i]){
      redmesh_cov_con <- secr::raster(mesh[[i]], cov_con[i])
      redmesh_cov_con <- raster::aggregate(redmesh_cov_con, fact = 3, fun = mean)
      con_name <- quo_name(cov_con[i])
      redmesh_cov_con_df <- as.data.frame(redmesh_cov_con) %>% dplyr::rename(!!con_name := layer)
      red_df_con <- cbind(as.data.frame(coordinates(redmesh_cov_con)), redmesh_cov_con_df, D = 1) %>% 
        filter_at(cov_con[i], all_vars(!is.na(.)))
      red_df <- red_df %>% left_join(red_df_con, by = c(x, y))
    }
    
    red_mesh[[i]] <- read.mask(data = red_df)
  }
  
  mesh <- red_mesh
  
  return(mesh)
}