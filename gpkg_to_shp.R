# function to convert a gpkg to a shp

gpkg_to_shp_to_zip <- function(path, file, type = "shp"){
  require(sf)
  require(tidyverse)
  input <- st_read(file.path(path, file))
  st_write(input, str_replace(file.path(path, file), "gpkg", type), delete_dsn=TRUE)
  pattern = str_split(file,"_to",simplify = T)[1,1]
  datestr = str_extract(file, "\\d{7}") %>% as.Date("%Y%j")
  datepat = paste("to",
                  lubridate::month(datestr, label=TRUE, abbr=FALSE),
                  lubridate::year(datestr),
                  sep="_")
  setwd(path)
  zfiles <- list.files(path, pattern = pattern)
  zip(zipfile= paste0(path,"/", pattern,"_", datepat, "_gpkg_shp.zip" ),
      files = zfiles)
}

gpkg_to_shp_to_zip(path="/home/a/Desktop/FIRED_stuff",
            file="fired_south_sudan_to2021182_events.gpkg")


# the_path <- "/home/a/Desktop/FIRED_stuff/upload_to_drive"
# to_update <- list.files(the_path, 
#                         pattern = "events.gpkg")
# 
# for(f in to_update) {
#   gpkg_to_shp_to_zip(path=the_path,
#                      file=f)
# }
