macro_design <- function(path_to_study_area, path_to_surveyed_pts = NULL, buffer = 10000,
                         lat_varname = "lat", lng_varname = "lng",
                         strat_var, strat_breaks = 3, strat_proportions = NULL,
                         initial_seed = 1234, n_existing_pts = 0, n_new_pts,
                         basemap = "Esri.WorldStreetMap"){

  ## get study area shape file
  nc_sf <- st_read(path_to_study_area)

  ## make strata
  nc_sf <- nc_sf %>% mutate(stratum_id = raster::cut(eval(rlang::sym(strat_var)), breaks = strat_breaks))

  ## get already surveyed areas
  already_surveyed <- load_already_surveyed_areas(path_to_surveyed_pts, buffer = buffer, lat_varname = lat_varname, lng_varname = lng_varname)
  
  ## set stratum proportions to area-proportional if not set manually
  if(is.null(strat_proportions)){
    strat_proportions <- nc_sf %>% st_set_geometry(NULL) %>%
      count(stratum_id) %>%
      add_tally(n, name = "nn") %>%
      mutate(props = n / nn) %>%
      arrange(stratum_id) %>% select(props) %>% unlist()
  }

  # set colours for stratification variable
  pal <- colorFactor("YlOrRd", domain = nc_sf$stratum_id)

  # generate lots of Halton points

  # these points will be randomly generated across the study region; we need enough to make sure
  # there are *at least* the desired number of survey points in each stratum. There is no harm in
  # generating many more points that you need, because any contiguous subsequence is also spatially
  # well distributed, and the sequence is ordered (earlier elements surveyed first)

  draw <- 10000
  set.seed(initial_seed)
  my_seeds <- round(c(runif(1, 1, 1e6), runif(1, 1, 1e6)))
  pts <- generate_Halton_pts(n = draw, seeds = my_seeds, bases = c(2,3))

  # scale and shift Halton points to fit the bounding box of the study region
  bb <- st_bbox(nc_sf)
  scale <- c(bb$xmax, bb$ymax) - c(bb$xmin, bb$ymin)
  shift <- c(bb$xmin, bb$ymin)

  pts <- t(t(pts) * scale + shift)
  pts <- data.frame(pts)
  colnames(pts) <- c("X", "Y")

  # make the data frame an sf object
  pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = st_crs(nc_sf)$proj4string)

  # remove any points that fall outside the study area
  pts_sf_in_studyarea <- lapply(st_intersects(pts_sf, nc_sf), length) %>% unlist()
  pts_sf <- pts_sf[pts_sf_in_studyarea == 1, ]

  # remove any points that are too close to already surveyed points
  if(!is.null(path_to_surveyed_pts)){
    pts_sf_in_already_surveyed_area <- lapply(st_intersects(pts_sf, already_surveyed$utm), length) %>% unlist()
    pts_sf <- pts_sf[pts_sf_in_already_surveyed_area == 0, ]
  }

  # work out what stratum each point falls into, and order the points sequentially both overall and within each stratum
  # (within-stratum order not necessary; global order used to determine survey ordering)
  pts_sf <- pts_sf %>% mutate(pt_id = 1:nrow(pts_sf),
                              stratum_id = nc_sf$stratum_id[st_intersects(pts_sf, nc_sf) %>% unlist()],
                              gridcell_id = nc_sf$Id[st_intersects(pts_sf, nc_sf) %>% unlist()]) %>%
    group_by(stratum_id) %>%
    mutate(ws_pt_id = row_number()) %>%
    ungroup() %>% cbind(st_coordinates(pts_sf))

  # add lat-long
  pts_sf_ll <- st_transform(pts_sf, crs = 4326)
  pts_sf_ll$lng <- st_coordinates(pts_sf_ll)[,1]
  pts_sf_ll$lat <- st_coordinates(pts_sf_ll)[,2]

  # add a variable indicating absolute order of adding points, this respects both Halton ordering
  # in pt_id and the desired number of points in each stratum
  obs_n_per_strat <- rep(0, length(strat_proportions))
  strat_ids <- c()
  for(i in 1:nrow(pts_sf_ll)){
    exp_n_per_strat <- i * strat_proportions
    this_strat <- which.min(obs_n_per_strat - exp_n_per_strat)
    strat_ids <- c(strat_ids, this_strat)
    obs_n_per_strat[this_strat] <- obs_n_per_strat[this_strat] + 1
  }
  order_df <- data.frame(ordered_id = 1:nrow(pts_sf_ll),
                         stratum_id = factor(levels(pts_sf_ll$stratum_id)[strat_ids],
                                             levels = levels(pts_sf_ll$stratum_id))) %>%
    group_by(stratum_id) %>%
    mutate(ws_pt_id = row_number()) %>% ungroup()

  pts_sf_ll <- left_join(pts_sf_ll, order_df, by = c("stratum_id", "ws_pt_id")) %>% arrange(ordered_id)

  # want to keep track of:
  # (a) total number of points generated, (b) the last N points added,
  # to plot (a) as existing points and (b) as the latest set of points.

  all_sampled_pts <- pts_sf_ll %>% filter(ordered_id <= (n_existing_pts + n_new_pts))

  existing_pts <- pts_sf_ll %>%
    filter(ordered_id <= n_existing_pts) %>%
    st_as_sf(coords = c("lng", "lat"), crs = 4326, remove = FALSE)

  new_pts <- anti_join(all_sampled_pts, existing_pts %>% st_set_geometry(NULL), by = "pt_id")

  ## plot

  m <- leaflet() %>%
    addTiles() %>%
    addProviderTiles(basemap) %>%
    addPolygons(data = st_transform(nc_sf, crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal(stratum_id)) %>%
    addScaleBar(position = "bottomleft") %>%
    addLegend(layerId = "myleg", colors = pal(levels(nc_sf$stratum_id)), opacity = 0.7,
              position = "bottomright", labels = paste0(levels(nc_sf$stratum_id), " (", 100*round(strat_proportions, 2), "%)")) %>%
    # add inset map
    addMiniMap(
      tiles = providers$Esri.WorldStreetMap,
      position = 'bottomleft',
      width = 150, height = 150,
      zoomLevelOffset = -4,
      toggleDisplay = TRUE) %>%
    addCircleMarkers(data = new_pts, layerId = ~pt_id, lng = ~lng, lat = ~lat, radius = 6, color = "red",
                     label = ~paste0("svy ", pt_id, " (lng ", round(lng,2), ", lat ", round(lat,2), ")")) %>%
    addCircleMarkers(data = existing_pts, layerId = ~pt_id, radius = 3, lng = ~lng, lat = ~lat, color = "blue",
                     label = ~paste0("svy ", pt_id, " (lng ", round(lng,2), ", lat ", round(lat,2), ")"))

  if(!is.null(path_to_surveyed_pts)){
    m <- m %>% addPolygons(data = already_surveyed$latlong, weight = 0, color = "black")
  }

  # measure of spatial balance
  N <- nrow(pts_sf_ll)
  n <- nrow(existing_pts) + nrow(new_pts)
  p <- rep(n/N, N)
  X <- pts_sf_ll %>% st_set_geometry(NULL) %>% select(X,Y) %>% as.matrix()
  s_exi <- existing_pts %>% st_set_geometry(NULL) %>% select(pt_id) %>% unlist()
  s_new <- new_pts %>% st_set_geometry(NULL) %>% select(pt_id) %>% unlist()
  s <- sort(c(s_exi, s_new))
  sbal <- round(BalancedSampling::sb(p, X, s), 3)

  return(list(m = m, sbal = sbal))

}
