# shift firedpy polygons
library(sf)
library(stringr)

posthoc_fixes <- function(filename, 
                          fix_offset=TRUE, 
                          add_crs=TRUE,
                          modis_crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs",
                          resolution = 463.31271653,
                          output_filename = str_replace(filename, ".gpkg", "_corrected.gpkg")){
  require(sf);require(stringr)
  d <- st_read(filename)
  
  if(fix_offset==TRUE) st_geometry(d) <- st_geometry(d)+ c(-resolution/2,0) + c(0, resolution/2)
  if(add_crs==TRUE) st_crs(d) <- st_crs(modis_crs)
  
  st_write(d, output_filename, delete_dsn = TRUE)
  return(d) # so you can test it
}

d<- posthoc_fixes(filename)
