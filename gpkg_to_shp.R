# function to convert a gpkg to a shp

gpkg_to_shp <- function(x, type = "shp"){
  require(sf)
  input <- st_read(x)
  st_write(input, str_replace(x, "gpkg", type))
}

gpkg_to_shp("/home/a/Desktop/FIRED_stuff/fired_ivory_coast_to2021182_events.gpkg")
