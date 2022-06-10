library(sf)
library(tidyverse)

world <- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  dplyr::select(NAME_EN, CONTINENT) %>%
  mutate(NAME_EN = str_to_lower(NAME_EN)%>% str_replace_all(" ","_"))

completed_countries <-c(filter(world, CONTINENT == "North America") %>% pull(NAME_EN),
                        filter(world, CONTINENT == "Europe") %>% pull(NAME_EN),
                        "uzbekistan", "afghanistan", "pakistan", "brazil", "zimbabwe",
                        "armenia", "georgia", "azerbaijan","mongolia", "eritrea",
                        "mali", "mauritania","guinea", "guinea-bissau", "sierra_leone",
                        "libya", "algeria", "tunisia","madagascar", "liberia",
                        "south_africa", "lesotho", "myanmar","laos", "cambodia",
                        "ivory_coast", "djibouti", "ethiopia", "somalia", "somaliland",
                        "ghana","togo", "benin", "niger","burkina_faso", "nigeria",
                        "colombia", "ecuador", "peru", "venezuela", "sudan","rwanda",
                        "south_sudan","central_african_republic", "namibia",
                        "cameroon", "gabon", "chad", "japan", "north_korea",
                        "uganda", "kenya", "burundi","australia", "malawi", "eswatini",
                        "malaysia", "brunei", "vietnam", "thailand", "indonesia",
                        "south_korea", "taiwan", "philippines", "papua_new_guinea",
                        "east_timor", "new_zealand", "nepal", "bhutan", "bangladesh",
                        "equatorial_guinea", "cameroon", "western_sahara", "sri_lanka",
                        "chile", "argentina", "uruguay", "paraguay","suriname", 
                        "bolivia", "guyana", "senegal", "morocco", "the_gambia","french guyana",
                        "turkmenistan", "people's_republic_of_china", "russia", "saudi_arabia", "jordan",
                        "lebanon", "iraq", "syria", "turkey", "oman", "yemen","tajikistan",
                        "united_arab_emirates", "india", "palestine", "israel",
                        "egypt", "iran", "kazakhstan", "kyrgyzstan", "botswana",
                        "mozambique", "republic_of_the_congo", "tanzania",
                        "angola", "zambia", "democratic_republic_of_the_congo")

lut_completed <- rep("Complete",length(completed_countries))
names(lut_completed) <- completed_countries


wrld<-world %>%
  mutate(completed = lut_completed[as.character(NAME_EN)])%>%
  replace_na(list(completed = "Coming soon")) %>%
  mutate(completed = replace(completed, NAME_EN == "antarctica" | NAME_EN == "greenland",
                             "Insufficient Fire\nActivity"))%>%
ggplot()+
  geom_sf(aes(fill = completed),lwd=0.10) +
  scale_fill_manual(values = c("skyblue", "red","grey95"),na.value = "grey95")+
  theme_void()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Fire Perimeter Datasets Available"," ")+
  theme(legend.position = c(.05,.3),
        legend.title = element_blank(),
        plot.title = element_text(hjust=.5, size=20),
        legend.justification = c(0,0));wrld

ggsave(wrld, filename = "/home/a/projects/firedpy/map_fig.png",
       bg="white",width=7, height=5, dpi=600)

# australia management
# ge<- st_read("/home/a/projects/firedpy/ref/individual_countries/georgia.gpkg")
# aus<-st_read("/home/a/Desktop/FIRED_stuff/STE_2021_AUST_GDA2020.shp") %>%
#   st_transform(crs = st_crs(ge))%>%
#   st_simplify(preserveTopology = F, dTolerance = 10000)
# 
# ggplot(aus)+
#   geom_sf()
# 
# for(i in aus$STE_NAME21[1:7]){
# aus %>%
#   filter(STE_NAME21 == i) %>%
#   st_write(paste0("/home/a/projects/firedpy/ref/australia_states/",
#                   str_to_lower(i) %>% str_replace_all(" ","_"), ".gpkg"),
#            delete_dsn=TRUE)
#   }
# 
# aus %>% filter(STE_NAME21 == "New South Wales" | STE_NAME21 == "Australian Capital Territory") %>%
#   summarise %>% 
#   st_write("/home/a/projects/firedpy/ref/australia_states/nsw_capital.gpkg", 
#            delete_dsn=TRUE)
# 
# st_read("/home/a/projects/firedpy/ref/west_australia.gpkg") %>% st_transform(crs=st_crs(ge)) %>%
#   summarise() %>%
#   st_simplify(preserveTopology = F, dTolerance = 10000)%>%
#   st_write("/home/a/projects/firedpy/ref/west_australia.gpkg", delete_dsn=T)


