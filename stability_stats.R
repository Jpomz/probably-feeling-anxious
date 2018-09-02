# Stability statistics
library(dplyr)
library(purrr)
library(ggplot2)
library(robustbase)
library(lmtest)
library(forcats)

dat <- readRDS("data/stability_results.RDS") %>%
  split(list(.$.id))

# heteroskedastic?
dat %>% map(~bptest(stab~pca1, data = .x))       
# all are

# linear or quadratic?
dat %>% map(~anova(lmrob(stab~pca1, data = .x),
                   lmrob(stab~pca1 +I(pca1^2), data = .x),
                   test = "Wald"))
# all quadratic
dat %>% map(~summary(lmrob(stab~pca1 +I(pca1^2), data = .x)))


fig_fun <- function(data, mainlab,
                    x = "pca1", y = "stab"){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    geom_point(aes(color = fct_reorder(Site, pca1)))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_discrete(h.start = 80) +
    labs(x = "Mining gradient",
         y = "stability metric",
         title = mainlab)+
    stat_smooth(method = "lmrob",
                formula = y~x+I(x^2),
                color = "black") 
}

fig_fun(data = dat$random, mainlab = "Random")
ggsave("figures/s~random.png")
fig_fun(data = dat$scaled, mainlab = "Scaled")
ggsave("figures/s~scaled.png")
fig_fun(data = dat$corr, mainlab = "Correlated")
ggsave("figures/s~correlated.png")
fig_fun(data = dat$s.c, mainlab = "Scaled & correlated")
ggsave("figures/s~s_c.png")
