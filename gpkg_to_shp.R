# function to convert a gpkg to a shp

gpkg_to_shp <- function(x, type = "shp"){
  require(sf)
  input <- st_read(x)
  st_write(input, str_replace(x, "gpkg", type), delete_dsn=TRUE)
  pattern = str_extract(x, "fired_[a-z]+")
  datepat = str_extract(x, "to\\d{7}")
  setwd("/home/a/Desktop/FIRED_stuff/")
  zfiles <- list.files("/home/a/Desktop/FIRED_stuff", pattern = pattern)
  zip(zipfile= paste0("/home/a/Desktop/FIRED_stuff/", pattern,"_", datepat, "_gpkg_shp.zip" ),
      files = zfiles)
}

gpkg_to_shp("/home/a/Desktop/FIRED_stuff/fired_benin_to2021182_events.gpkg")
