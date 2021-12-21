# counting reburns
library(terra)
library(tidyverse)
terra::rast("/home/a/Desktop/FIRED_stuff/h20v09.nc") ->nc



bands <- seq(from=3, to =250, by = 12)

result<- data.frame(year=NA, percent_reburn=NA, count_reburn = NA, count_burned = NA)
counter<-1
for(b in bands[1:20]){
print(counter)
year<-nc[[b:(b+11)]] 

year[year>0] <- 1 
year[year<0] <- 0

tab<- year %>%
  sum() %>%
  terra::freq() %>%
  as_tibble() %>%
  filter(value>0) 

burns <- tab %>% filter(value == 1) %>% pull(count)
reburns <- tab %>% filter(value > 1) %>% pull(count) %>% sum()

result[counter,1] <- 2000+counter
result[counter,2] <- 100*(reburns/sum(c(burns, reburns)))
result[counter,3] <- reburns
result[counter,4] <- sum(c(burns, reburns))
counter<- counter+1
terra::tmpFiles(remove=TRUE)
gc()

}

table <- result %>%
  mutate(reburn_ha = count_reburn * (463^2) * 0.0001,
         reburn_km2 = count_reburn * (463^2) / 1000000) %>%
  pivot_longer(cols = percent_reburn:reburn_km2, 
               names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(mean = mean(value),
            std = sd(value),
            min = range(value)[1],
            max = range(value)[2])

