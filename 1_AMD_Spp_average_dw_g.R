# 2_AMD_Spp_average_dw_g

library(plyr)
library(tidyverse)

bugs <- read_csv("data/raw data/AMD estimated invert dw g.csv")
fish <- read_csv("data/raw data/estimated fish dw g.csv")

# remove Burcoa site
bugs <- bugs %>%
  filter(site != "Burcoa")
# add area column
# one surber area = 0.24 * 0.25 = 0.06 m^2
bugs$area <- 0.06
# change Label to taxa for row binding below
bugs <- rename(bugs, taxa = Genus)
# drop variables
bugs <- bugs %>% select(taxa, site, surber, dw, area)

bugs <- bugs %>%
  group_by(site,  taxa) %>%
  summarize(avg.dw = mean(dw),
            count = n(),
            area = mean(area),
            density = (count / 3) * 16.66667)


fish <- fish %>% 
  select(taxa, site, dw, area)
fish <- fish %>%
  group_by(site, taxa) %>%
  summarize(avg.dw = mean(dw),
            count = n(),
            area = mean(area),
            density = count / area)


# bind rows
full.data <- bind_rows(bugs, fish)

# calculate average dw, density, area, per site and taxa
dataset <- full.data %>% 
  group_by(site) %>%
  mutate(tot.ab = sum(density),
         rel.ab = density / tot.ab) %>%
  ungroup()

saveRDS(dataset, "data/AMD fish invert dw abundance.RDS")
