library(tidyverse)
library(sf)
modis_tiles <- st_read("/home/a/data/background/modis_grid/modis_sinusoidal_grid_world.shp")

wb<- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  st_transform(crs = st_crs(modis_tiles))

modtiles <- st_read("/home/a/projects/firedpy/ref/modis_grid.gpkg")
st_crs(modtiles) <- st_crs(modis_tiles)

p1<- ggplot() +
  geom_sf(data=wb, fill="burlywood", lwd=0.25) +
  geom_sf(data =modtiles, fill="transparent",color="black", lwd=0.25) +
  theme(panel.background = element_rect(fill = "transparent", color="black"));p1

ggsave(p1, filename="/home/a/Desktop/modis_tile_plot.png", 
       dpi=300,
       width = 7.5, 
       height=4)

xx<-st_read("/home/a/projects/firedpy/fired_italy_to2021182_events.gpkg")
yy<-st_read("/home/a/projects/firedpy/fired_italy_to2021182_daily.gpkg")


p2<-bind_rows(
  yy %>% filter(id == 192) %>%
    dplyr::select(id, Date=date) %>% 
    group_by(Date) %>%
    summarise() %>%
    mutate(product = "A. Daily"),
  xx %>% filter(id == 192) %>% dplyr::select(id, Date=ig_date)%>%  
  mutate(product = "B. Event"))%>%
  ggplot() +
  geom_sf(aes(fill = Date), lwd=0.25) +
  theme_minimal()+
  scale_fill_viridis_d()+
  facet_wrap(~product) +
  theme(legend.position = "left",
        axis.text = element_blank(),
        plot.background = element_rect(color="black"),
        panel.background = element_rect(color="black"))

ggsave(p2, filename="/home/a/projects/firedpy/daily_example.png", dpi=300, height=3.5, width=4.5)
