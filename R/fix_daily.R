# fix daily polygons not dissolving
library(sf)
library(tidyverse)
# devtools::install_github("hadley/multidplyr")
library(multidplyr)
dfile <-"/home/a/data/fire/fired/wa_or_ca/fired_wa_or_ca_to2021305_daily.gpkg"

dsf <- st_read(dfile) %>%
  mutate(date = as.Date(date)) %>%
  filter(date > as.Date("2021-07-01"))

# Creating 4-core cluster
cl <- new_cluster(4)
# cl <- default_cluster()
cluster_library(cl, "dplyr")

library(doParallel)

library(foreach)
corz<-4
registerDoParallel(corz)



fix_daily <- function(dsf, cl){

    att<- dsf %>%
      st_set_geometry(NULL) %>%
      group_by(id, date) %>%
      summarise_all(first) %>%
      dplyr::select(-did, -x, -y)
    
    ids<-unique(dsf$id)
    ns <- seq(1, length(ids), by=100)
    
    geom <- foreach(i = ns, .combine = bind_rows)%dopar%{
      system(paste("echo", i, "/", length(ids)))
      
      idsub <- ids[i:(i+99)]
      
      dsf %>%
        filter(id %in% ids) %>%
        dplyr::select(id, date) %>%   
        group_by(id, date) %>%
        summarise() %>%
        ungroup()
                  }
    
    left_join(geom,att, by = c("id", "date"))->x
    return(x)
}

st_write(x,"/home/a/data/fire/fired/wa_or_ca/fired_wa_or_ca_to2021305_daily_fixed.gpkg" )
