library(sf)
library(tidyverse)

world <- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  dplyr::select(NAME_EN, CONTINENT) %>%
  mutate(NAME_EN = str_to_lower(NAME_EN)%>% str_replace_all(" ","_"))


plot()
completed_countries <-c(filter(world, CONTINENT == "North America") %>% pull(NAME_EN),
                        filter(world, CONTINENT == "Europe") %>% pull(NAME_EN),
                        "uzbekistan", "afghanistan", "pakistan", "brazil",
                        "chile", "argentina", "uruguay", "paraguay","suriname", 
                        "bolivia", "guyana", "senegal", "morocco", "the_gambia","french guyana",
                        "turkmenistan", "people's_republic_of_china", "russia", "saudi_arabia", "jordan",
                        "lebanon", "iraq", "syria", "turkey", "oman", "yemen","tajikistan",
                        "united_arab_emirates", "india", "palestine", "israel",
                        "egypt", "iran", "kazakhstan", "kyrgyzstan")
lut_completed <- rep("complete",length(completed_countries))
names(lut_completed) <- completed_countries


wrld<-world %>%
  mutate(completed = lut_completed[as.character(NAME_EN)])%>%
ggplot()+
  geom_sf(aes(fill = completed),lwd=0.25) +
  scale_fill_manual(values = c("red", "white"),na.value = "grey90")+
  theme_void()+
  ggtitle("Completed Countries")+
  theme(legend.position = c(.1,.3),
          plot.title = element_text(hjust=.5),
        legend.justification = c(0,0));wrld

ggsave(wrld, filename = "/home/a/projects/firedpy/completed_countries_plot.png",
       bg="white",width=7, height=5, dpi=300)
# 
# Coterminous USA + Alaska
# US plus Canada
# Hawaii
# All the countries in the Carribean
# Mexico and Central America
# 
# South America
# 
# Bolivia
# Argentina
# Northern South America (Suriname, French Guiana, Guyana)
# Chile
# Uruguay
# Brazil
# 
# 
# Europe (November 2000 to July 2021)
# 
# Northern Europe (ICELAND, SWEDEN, NORWAY, and DENMARK)
# Russia
# Italy
# Spain & Portugal
# Western Europe (FRANCE, GERMANY, POLAND, SWITZERLAND, BELGIUM, NETHERLANDS, LUXEMBOURG and AUSTRIA)
# Central to Southern Europe (ESTONIA, LATVIA, LITHUANIA, BELARUS, UKRAINE, CZECH REPUBLIC, SLOVAKIA, HUNGARY, ROMANIA, BULGARIA, MONTENEGRO, BOSNIA, TURKEY, REPUBLIC OF MOLDOVA, SERBIA, ALBANIA, SLOVENIA, and NORTH MACEDONIA)
# Greece
# UK and Ireland
# 
# Africa
# 
# Senegal
# Gambia
# Morocco
# 
# Asia
# 
# China
# India
# Central Asia (Turkmenistan, Kazakhstan, Uzbekistan, Kyrgystan, Tajikistan, Afghanistan, and Pakistan)
# Middle East (Saudi Arabia, Quatar, Oman, Yemen, United Arab Emirates, Iraq, Jordan, Syria, Israel, Palestine, Lebanon, Egypt)
