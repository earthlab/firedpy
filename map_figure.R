library(sf)
library(tidyverse)

world <- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  dplyr::select(NAME_EN, CONTINENT) %>%
  mutate(NAME_EN = str_to_lower(NAME_EN)%>% str_replace_all(" ","_"))


plot()