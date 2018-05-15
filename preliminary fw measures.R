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



#library(Hmisc)

fwpointrange <- function(data, x = "pca1", y, ylab = NULL){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x),
                 fun.ymax = function(x) mean(x) + sd(x),
                 geom = "pointrange") +
    theme_bw() +
    theme(axis.title = element_text(size = 20)) +
    if (is.null(ylab))
      labs(x = "Mining gradient", y = y) 
  else
    labs(x = "Mining gradient", y = ylab) 
}

fwpointrange(fw.measures, y = "L") +
  geom_smooth(method = "lm")
ggsave("figures/L.png")
fwpointrange(fw.measures, y = "C") +
  geom_smooth(method = "lm")
ggsave("figures/C.png")
fwpointrange(fw.measures, y = "Gensd", ylab = "SD(Gen)") +
  geom_smooth(method = "lm")#, formula = y ~ x + I(x^2))
ggsave("figures/Gsd.png")
fwpointrange(fw.measures, y = "Vulsd", ylab = "SD(Vul)") +
  geom_smooth(method = "lm")#, formula = y ~ x + I(x^2))
ggsave("figures/Vsd.png")
fwpointrange(fw.measures, y = "B") +
  geom_smooth(method = "lm")
fwpointrange(fw.measures, y = "I") +
  geom_smooth(method = "lm")
fwpointrange(fw.measures, y = "T") +
  geom_smooth(method = "lm")
fwpointrange(fw.measures, y = "Maxsim") +
  geom_smooth(method = "lm")
fwpointrange(fw.measures, y = "max.TL") +
  geom_smooth(method = "lm")

test <- AIC(lm(C ~ pca1, fw.measures), lm(C~ pca1 + I(pca1^2), fw.measures))
    

summary(lm(C~ pca1 + I(pca1^2), fw.measures))
summary(lm(C~ pca1 , fw.measures))
