# Preliminary food web measures

library(plyr)
library(tidyverse)
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

# network data
dat <- readRDS("data/500 bernouli trials.rds")


summary(dat)
fw.measures <- llply(dat, function (x){
  map(x, Get.web.stats)
})

fw.measures <- llply(fw.measures, function (x){
  ldply(x)
}) %>% ldply()

pca.axis <- readRDS("data/pca axis 26 sites.rds")[,c(4, 2)]
fw.measures <- left_join(fw.measures, pca.axis, by = c(".id" = "site"))

fwline <- function(data, x = "pca1", y, ylab = NULL){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    geom_point() +
    theme_classic() +
    stat_smooth(method = "lm") +
    if (is.null(ylab))
      labs(x = "PC1", y = y)
  else
    labs(x = "PC1", y = ylab)
}

fwline(fw.measures, y = "S")
fwline(fw.measures, y = "L")
fwline(fw.measures, y = "C")
fwline(fw.measures, y = "B")
fwline(fw.measures, y = "I")
fwline(fw.measures, y = "T")
fwline(fw.measures, y = "Gensd")
fwline(fw.measures, y = "Vulsd")
fwline(fw.measures, y = "max.TL")

fw.measures %>% group_by(.id, pca1) %>%
  summarize(L.m = mean(L), L.sd = sd(L),
            C.m = mean(C), C.sd = sd(C),
            G.m = mean(Gensd), G.sd = sd(Gensd, na.rm = TRUE),
            V.m = mean(Vulsd), V.sd = sd(Vulsd, na.rm = TRUE)) %>%
  arrange(pca1) %>% View








