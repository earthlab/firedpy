# fix daily polygons not dissolving
library(sf)
library(tidyverse)
# devtools::install_github("hadley/multidplyr")
# library(multidplyr)
library(doParallel)
library(foreach)

# functions =======================
fix_daily <- function(dsf, cl=4){

    att<- dsf %>%
      st_set_geometry(NULL) %>%
      group_by(id, date) %>%
      summarise_all(first) %>%
      dplyr::select(-did, -x, -y)

    ids<-unique(dsf$id)
    ns <- seq(1, length(ids), by=10)

    geom <- foreach(i = ns, .combine = bind_rows)%dopar%{
      system(paste("echo", i, "/", length(ids)))

      idsub <- ids[i:(i+9)]

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

corz<-detectCores()-1
registerDoParallel(corz)

dfile <-"fired_uscan_to2021121_daily.gpkg"
dfile_out <- paste0(str_split(dfile, "\\.")[[1]][1], "_fixed.gpkg")


ids <- st_read(dfile, query="SELECT id FROM fired_uscan_to2021121_daily", quiet=TRUE) %>%
  pull(id) %>%
  unique() %>%
  sort()

gc()

iterators <- seq(1, length(ids), by = 100)

dir.create("tmp")
# result <- list()
for(Z in 1:length(iterators)){
  tmpfl <- paste0("tmp/tmp",Z,".gpkg")
  print(paste("chunk", Z, "of", length(iterators)))

  if(!file.exists(tmpfl)){
    
    if(Z != length(iterators)){
      id_subset <- ids[iterators[Z]:(iterators[Z+1]-1)]
    }else{id_subset <- ids[iterators[Z]:length(ids)]}
    
    query <- paste("SELECT * FROM fired_uscan_to2021121_daily WHERE id >=",
                   min(id_subset),
                   "AND id <=", max(id_subset))
  
    dsf <- st_read(dfile, query=query) %>%
      mutate(date = as.Date(date))
  
    fix_daily(dsf, cl) %>%
      st_write(tmpfl, delete_dsn=T)
    gc()
    system(paste0("aws s3 cp ", tmpfl, " s3://earthlab-amahood/", tmpfl))
  }
}

fl <- list.files("tmp", full.names=T)

lapply(fl, st_read) %>%
  bind_rows() %>%
  st_write(dfile_out)
system(paste0("aws s3 cp ", dfile_out, " s3://earthlab-amahood/", dfile_out))
