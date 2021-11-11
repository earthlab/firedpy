# function to convert a gpkg to a shp

gpkg_to_shp_to_zip <- function(path, file, type = "shp"){
  require(sf)
  input <- st_read(file.path(path, file))
  st_write(input, str_replace(x, "gpkg", type), delete_dsn=TRUE)
  pattern = str_extract(file, "fired_[a-z]+")
  datepat = str_extract(file, "to\\d{7}")
  setwd(path)
  zfiles <- list.files(path, pattern = pattern)
  zip(zipfile= paste0(path,"/", pattern,"_", datepat, "_gpkg_shp.zip" ),
      files = zfiles)
}

gpkg_to_shp_to_zip(path="/home/a/Desktop/FIRED_stuff",
            file="fired_victoria_to2021182_events.gpkg")
