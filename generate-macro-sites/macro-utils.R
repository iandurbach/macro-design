load_already_surveyed_areas <- function(filename, buffer, lat_varname, lng_varname){

  if(!is.null(filename)){
    x <- read_csv(filename)
    
    # rename lat_varname and lng_varname to lat and lng
    x <- x %>% dplyr::rename(lat = !!lat_varname, 
                      lng = !!lng_varname)
    
    # remove any missing values on lat/long
    x <- x %>% filter(!is.na(lng)) %>% filter(!is.na(lat))
         
    # turn into sf object
    already_sf <- st_as_sf(x, coords = c("lng", "lat"), crs = 4326)
    # project for buffering
    already_sf_utm <- st_transform(already_sf, crs = st_crs(nc_sf)$proj4string)
    # add 10km buffer
    already_sf_buffered <- st_buffer(already_sf_utm, dist = buffer)
    # convert back to ll
    already_sf_ll <- st_transform(already_sf_buffered, crs = 4326)
  }else{
    already_sf <- NA
    already_sf_utm <- NA
    already_sf_ll <- NA
    already_sf_buffered <- NA
  }

  return(list(cams_ll = already_sf, cams_utm = already_sf_utm,
              latlong = already_sf_ll, utm = already_sf_buffered))

}
