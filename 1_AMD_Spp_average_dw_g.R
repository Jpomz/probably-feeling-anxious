# 1_AMD_Spp_average_dw_g

library(plyr)
library(tidyverse)

bugs <- read_csv("data/raw data/AMD estimated invert dw g.csv")
fish <- read_csv("data/raw data/estimated fish dw g.csv")

# add rows for total estimated fish calculated in 1a_k_pass_removal.R script 
# Only kiwi site had additional fish
# --> add 2 anguilla, 3 gobiomorphus, 2 galaxias maculatus

# calculate mean dw values for additional fish
mean.k.ang <- fish %>%
  filter(site =="Kiwi", taxa =="Anguilla.australis") %>%
  summarize(dw = mean(dw))

mean.k.gob <- fish %>%
  filter(site =="Kiwi", taxa =="Gobiomorphus.huttoni") %>%
  summarize(dw = mean(dw))

mean.k.mac <- fish %>%
  filter(site =="Kiwi", taxa =="Galaxias.maculatus") %>%
  summarize(dw = mean(dw))

# make data frame of extra observations
plus.fish <- data.frame(site = "Kiwi",
                       date_coll = "2.3.16",
                       taxa = c(rep("Anguilla.australis",2),
                                rep("Gobiomorphus.huttoni", 3),
                                rep("Galaxias.maculatus",2)),
                       reach_width = 2.5, 
                       reach_length = 20,
                       ln_a = NA, 
                       b = NA, 
                       base = NA,
                       log = NA,
                       dw = as.numeric(c(rep(mean.k.ang, 2),
                              rep(mean.k.gob, 3),
                              rep(mean.k.mac, 2))),
                       area = 50,
                       stringsAsFactors = FALSE)

fish <- bind_rows(fish, plus.fish)

# cleanup bugs data
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

# calculate average dw, density, and area by site and taxa
dataset <- full.data %>% 
  group_by(site) %>%
  mutate(tot.ab = sum(density),
         rel.ab = density / tot.ab) %>%
  ungroup()

saveRDS(dataset, "data/AMD_fish_invert_dw_abundance.RDS")
