# figures for data paper

# setup ========================================================================

library(tidyverse)
library(sf)
library(ggpubr)
library(ggsn)
library(raster)
library(egg)

data_dir <- "/home/a/projects/firedpy/data_paper_files"
template <- raster::raster(file.path(data_dir, "template.tif"))

# rim fire comparison figure ===================================================

FIRED_p <- st_read(file.path(data_dir, "rim_fired.gpkg")) %>%
  st_set_crs(raster::crs(template, asText=TRUE)) %>%
  mutate(product = "b. FIRED (3 Events)") %>%
  dplyr::select(id, product)
gfa_p <- st_read(file.path(data_dir, "rim_gfa.gpkg")) %>%
  st_transform(crs = crs(FIRED_p)) %>%
  mutate(product = "c. GFA (14 Events)") %>%
  dplyr::select(id = fire_ID, product)
mtbs_p <- st_read(file.path(data_dir, "rim_mtbs.gpkg")) %>%
  st_transform(crs = crs(FIRED_p)) %>%
  mutate(product = "a. MTBS (1 Event)",
         id=1) %>%
  dplyr::select(id, product)
gwis_p <- list.files(file.path(data_dir, "gwis_1_8-9-10_2013_rim_fire"), full.names = TRUE) %>%
  lapply(st_read) %>%
  bind_rows() %>%
  st_transform(crs = crs(FIRED_p)) %>%
  mutate(product = "d. GWIS (48 Events)")  %>%
  filter(Type == "FinalArea") %>%
  dplyr::select(id=Id, product) 

all4 <- bind_rows(FIRED_p, gfa_p) %>%
  bind_rows(mtbs_p) %>% bind_rows(gwis_p)

event_figure <- ggplot(all4) +
  geom_sf(aes(fill=as.factor(id)), color="transparent") +
  facet_wrap(~product, nrow=2) +
  theme_classic() +
  scale_fill_identity()+
  theme(legend.position = "none",
        strip.text = element_text(size=15),
        panel.grid.major = element_line(color = "grey90"),
        panel.border =element_rect(size=1, fill="transparent"),
        axis.title = element_blank()) 

tag_pool <- all4 %>% arrange(product) %>% pull(product) %>% str_sub(4,99) %>% unique()

event_figure_t <- egg::tag_facet(event_figure, open="", close="",
                                       fontface = "bold", size=8)%>%
  egg::tag_facet(open="", close="", 
                 vjust = -.5, hjust = 1.1,
                 fontface = "plain", tag_pool = tag_pool,
                 x=Inf, y=-Inf, size=6)

ggsave(event_figure_t, filename = file.path(data_dir, "figure_1_rim_event_figure.png"), 
       bg="white", height = 8, width=12)
ggsave(event_figure_t, filename = file.path(data_dir, "figure_1_rim_event_figure.pdf"),
       bg="white", height = 8, width=12, dpi=600)

# table comparing rim fire to other products ===================================

rim_areas_main <- all4 %>%
  mutate(area = st_area(.) %>% units::drop_units(),
         area = area/1000000) %>%
  arrange(area %>% desc) %>%
  slice(1:4) %>%
  st_set_geometry(NULL)

rim_areas_secondary <- all4 %>%
  mutate(area = st_area(.) %>% units::drop_units(),
         area = area/1000000) %>%
  arrange(area %>% desc) %>%
  slice(5:97) %>%
  st_set_geometry(NULL) %>%
  group_by(product) %>%
  summarise(secondary_event_area = sum(area),
            n_secondary_events = n()) %>%
  ungroup() %>%
  bind_rows(data.frame(product = "A. MTBS (1 Event)", secondary_event_area = 0, n_secondary_events = 0))

left_join(rim_areas_main, rim_areas_secondary) %>%
  mutate(total_area = area + secondary_event_area) %>%
  mutate_if(is.numeric, round,1)%>%
  write_csv(file.path(data_dir, "mtbs_comparison.csv"))


# figure for small fires =======================================================

dates <- st_read(file.path(data_dir, "moonshine_sour_mtbs.gpkg")) %>% 
  st_set_geometry(NULL) %>%
  pull(Ig_Date)
fired_ms <- st_read(file.path(data_dir, "moonshine_sour_fired.gpkg")) %>%
  st_set_crs(raster::crs(template, asText=TRUE))%>%
  mutate(product = "A. FIRED 5 pixels, 11 days\n(1 Event)") %>%
  dplyr::select(id, product)

mtbs_ms <- st_read(file.path(data_dir, "moonshine_sour_mtbs.gpkg")) %>%
  st_transform(crs = crs(fired_ms)) %>%
  mutate(product = "A. MTBS\n(2 Events)",
         id=c(3,2)) %>%
  dplyr::select(id, product)

fired_t5s1 <- st_read(file.path(data_dir, "moonshine_sour_fired_t5s1.gpkg"))  %>%
  st_set_crs(raster::crs(template, asText=TRUE))%>%
  mutate(product = "B. FIRED 1 pixel, 5 days\n(1 Event)") %>%
  filter(ig_date > as.Date(dates[2]),
         ig_date < as.Date(dates[1])+5) %>%
  dplyr::select(id, product)

ints_fired <- st_intersects(fired_t5s1, mtbs_ms, sparse = F)
fired_s1 <- fired_t5s1[rowSums(ints_fired)>0,] %>%
  mutate(product = paste0("B. FIRED 1 pixel, 5 days\n(",nrow(.)," Events)"))
  


gwis_msr <- st_read(file.path(data_dir, 
                            "gwis_1_1-2-3-12_2007_1-2-3_2008_sour_moonshine", 
                            "gwis_1_3_2007_sour_moonshine.gpkg")) %>% 
  st_transform(crs = crs(fired_ms)) %>%
  mutate(product = "D. GWIS (x Events)")  %>%
  filter(Type == "FinalArea",
         IDate > as.Date(dates[2])) %>%
  dplyr::select(id=Id, product) 

ints_gwis <- st_intersects(gwis_msr, mtbs_ms, sparse = F)
gwis_ms <- gwis_msr[rowSums(ints_gwis)>0,] %>%
  mutate(product = paste0("D. GWIS\n(",nrow(.)," Events)"))



gfa_msr <- st_read(file.path(data_dir, "moonshine_sour_gfa.gpkg")) %>%
  st_transform(crs = crs(fired_ms)) %>%
  mutate(product = "C. GFA\n(14 Events)") %>%
  filter(start_date > as.Date(dates[2]),
         start_date < as.Date(dates[1])+5)%>%
  dplyr::select(id = fire_ID, product)

ints_gfa <- st_intersects(gfa_msr, mtbs_ms, sparse = F)
gfa_ms <- gfa_msr[rowSums(ints_gfa)>0,] %>%
  mutate(product = paste0("C. GFA\n(",nrow(.)," Events)"))

all4_ms <- bind_rows(fired_ms, gfa_ms) %>%
  bind_rows(fired_s1) %>% bind_rows(gwis_ms)

small_event_figure <- ggplot(all4_ms) +
  geom_sf(aes(fill=as.factor(id)), color="transparent") +
  geom_sf(data = mtbs_ms %>% dplyr::select(-product), color = "black", fill="transparent")+
  facet_wrap(~product, nrow=2) +
  theme_classic() +
  scale_fill_identity()+
  theme(legend.position = "none",
        strip.text = element_text(size=15),
        panel.grid.major = element_line(color = "grey90"),
        panel.border =element_rect(size=1, fill="transparent"),
        axis.title = element_blank())

tag_pool_ms <- all4_ms %>% arrange(product) %>% pull(product) %>% str_sub(4,99) %>% unique()

small_event_figure_t <- egg::tag_facet(small_event_figure, open="", close="",
                                 fontface = "bold", size=8)%>%
  egg::tag_facet(open="", close="", 
                 vjust = -.5, hjust = 1.1,
                 fontface = "plain", tag_pool = tag_pool_ms,
                 x=Inf, y=-Inf, size=6)

ggsave(small_event_figure_t, filename = file.path(data_dir, "figure_2_small_event_figure.png"), bg="white", height = 8, width=8)
ggsave(small_event_figure_t, filename = file.path(data_dir, "figure_2_small_event_figure.pdf"),
       bg="white", height = 8, width=8, dpi=600)

# how to st_intersects

# countries%>%
#   filter(st_intersects(car)[[1]])
# 
# ints <- st_intersects(gwis_ms, mtbs_ms, sparse = F)
# gwis_ms[rowSums(ints)>0,]

# Daily figure =================================================================

rim<-st_read(file.path(data_dir, "rim_fire.gpkg")) %>%
  st_set_crs(raster::crs(template, asText=TRUE))%>%
  mutate(source = "A. FIRED Event",
         last_date = lubridate::ymd(last_date)) %>%
  dplyr::select(date=last_date,source, geometry=geom)
  

rim_daily <- st_read(file.path(data_dir,"rim_daily.gpkg"), quiet=T) %>%
  mutate(date = lubridate::ymd(date))%>%
  st_set_crs(raster::crs(template, asText=TRUE))

filez <- list.files(file.path(data_dir,"rim_daily_ir"), pattern = ".shp$", full.names = TRUE)

list_of_polygonz <- lapply(filez, st_read) 

daily_ir <- bind_rows(list_of_polygonz) %>%
  dplyr::select(CaptDate) %>%
  dplyr::filter(CaptDate > as.Date("2012-01-02")) %>%
  st_transform(crs=st_crs(rim_daily))


both<- daily_ir %>% 
  arrange(desc(CaptDate)) %>%
  rename(date = CaptDate) %>%
  mutate(source = "C. Daily IR Flights") %>%
  bind_rows(rim_daily %>%
              rename(geometry = geom) %>%
              dplyr::select(date) %>% 
              mutate(source = "B. FIRED Daily", date = lubridate::ymd(date))) %>%
  bind_rows(rim) %>%
  dplyr::rename(Date=date)



panel<-ggplot(both, aes(fill = Date)) +
  geom_sf(color = "transparent")+
  scale_fill_date(low = "yellow", high = "blue")+
  theme_classic()+
  facet_wrap(~source) +
  theme(strip.text = element_text(size=15),
        panel.grid.major = element_line(color = "grey90"),
        panel.border =element_rect(size=1, fill="transparent"),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        axis.title = element_blank())

tag_pool_d <- both %>% arrange(source) %>% pull(source) %>% str_sub(4,99) %>% 
  unique() %>% str_replace(" ", "\n")

panel_t<- egg::tag_facet(panel, open="", close="",
                  fontface = "bold", size=8)%>%
  egg::tag_facet(open="", close="", 
                 vjust = -.3, hjust = 1.1,
                 fontface = "plain", tag_pool = tag_pool_d,
                 x=Inf, y=-Inf, size=6)
ggsave(panel_t, filename=file.path(data_dir, "figure_4_daily_illustration.png"), bg="white", height=5, width=10)
ggsave(panel_t, filename=file.path(data_dir, "figure_4_daily_illustration.pdf"), 
       bg="white", height=5, width=10, dpi = 600)


# Map Figure ===================================================================
modis_grid <- st_read("/home/a/projects/firedpy/ref/modis_grid.gpkg") %>%
  st_set_crs(value =st_crs(FIRED_p))
world <- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  st_transform(crs=st_crs(FIRED_p))

gridplot <- ggplot() +
  geom_sf(data = world, fill="burlywood", lwd=0.25) +
  geom_sf(data = modis_grid, fill="transparent", lwd=0.25) +
  theme_classic() +
  theme(plot.background = element_rect(color="black"))

ggsave(gridplot, filename = file.path(data_dir, "figure_3_modis_grid.png"),
       bg="white", height = 4, width=7.5, dpi=600)
ggsave(gridplot, filename = file.path(data_dir, "figure_3_modis_grid.pdf"),
       bg="white", height = 4, width=7.5, dpi=600)
