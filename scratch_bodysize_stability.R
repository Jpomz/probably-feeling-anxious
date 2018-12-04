# body size ~ stability?
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
stab <- readRDS("data/stability_results.RDS") %>%
  split(list(.$.id))
dw <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
gradient <- readRDS("data/pca_axis_26_sites.RDS")

dw.summary <- dw %>%
  group_by(site) %>%
  summarize(min.dw = log10(min(avg.dw)),
            mean.dw = log10(mean(avg.dw)),
            median.dw = log10(median(avg.dw)),
            max.dw = log10(max(avg.dw)),
            quant.05 = log10(quantile(avg.dw,
                                      probs = 0.05)),
            quant.95 = log10(quantile(avg.dw,
                                      probs = 0.95))) %>%
  left_join(gradient[,c(4,2)], by = "site")

dw.summary <- dw.summary %>%
  gather("variable", "value", 2:7)

  
ggplot(dw.summary, aes(y = value, x = pca1, color = variable)) +
  geom_point() +
  theme_bw() +
  stat_smooth(se = FALSE)

z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dw.summary <- dw.summary %>%
  group_by(variable) %>%
  mutate(z.score = z_score(value))

ggplot(dw.summary, aes(y = z.score, x = pca1,
                       color = variable)) +
  geom_point() +
  theme_bw() +
  stat_smooth(se = FALSE)

dw.summary %>%
  filter(pca1 > -1) %>%
ggplot(aes(y = z.score, x = pca1, color = variable)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  #stat_smooth(se = FALSE) +
  NULL

ggplot(stab$corr, aes(y = stab, x = pca1)) +
  geom_point(aes(color = Site)) +
  theme(legend.position = "none")

dw.summary %>%
  filter(variable == "quant.95",
         value < -2) %>%
  left_join(stab$corr, by = "pca1") %>%
  ggplot(aes(y = stab, x = value, color = pca1)) +
  geom_point() +
  theme_bw()+
  theme(legend.position = "none") +
  stat_summary(fun.y = quantile,
               fun.args = list(probs = c(0.05, 0.5, 0.95)),
               geom = "point",
               shape = 3,
               size = 5, 
               color = "red") +
  stat_quantile()
