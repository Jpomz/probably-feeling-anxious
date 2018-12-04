# Stability statistics
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(robustbase)
library(lmtest)
library(forcats)
library(gridExtra)

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


stab_fig_fun <- function(data, mainlab,
                    x = "pca1", y = "stab"){
  ggplot(data, aes(x = -data[[x]], y = data[[y]])) +
    geom_point(aes(color = fct_reorder(Site, pca1)))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_discrete(h.start = 80) +
    labs(x = "AMD gradient",
         y = "stability metric",
         title = mainlab)+
    stat_smooth(method = "lmrob",
                formula = y~x+I(x^2),
                color = "black") 
}

# Fig 3A-D ####
fig.3a <- stab_fig_fun(data = dat$random,
                       mainlab = "Random") +
  annotate(geom = "text", x = -5.6, y = 0.22, 
           label = "A", size = 6)
fig.3b <- stab_fig_fun(data = dat$scaled,
                       mainlab = "Scaled") +
  annotate(geom = "text", x = -5.6, y = 0.075,
           label = "B", size = 6)
fig.3c <- stab_fig_fun(data = dat$corr,
                       mainlab = "Correlated") +
  annotate(geom = "text", x = -5.6, y = 0.17,
           label = "C", size = 6)
fig.3d <- stab_fig_fun(data = dat$s.c,
                       mainlab = "Scaled & correlated") +
  annotate(geom = "text", x = -5.6, y = 0.11,
           label = "D", size = 6)

png(filename = "figures/fig_3.png",
     width = 190, height = 190, units = "mm", res =300)
grid.arrange(fig.3a, fig.3b, fig.3c, fig.3d, ncol = 2)
dev.off()


# # exploring other figs ####
# dat.df <- ldply(dat)
# 
# ggplot(dat.df, aes(y = stab, x = pca1, color = .id)) +
#   geom_point() +
#   stat_smooth(method = "lmrob",
#               formula = y~x+I(x^2))
# 
# dat.df %>% filter(.id == "random" | .id == "corr") %>%
# ggplot(aes(y = stab, x = pca1, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
# 
# # random vs sc
# dat.df %>% filter(.id == "random" | .id == "s.c") %>%
#   ggplot(aes(y = stab, x = pca1, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
# 
# # only look at unimpacted sites
# dat.df %>% filter(pca1 < -1,
#   .id == "random" | .id == "s.c") %>%
#   ggplot(aes(y = stab, x = pca1, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
# 
# # stability ~ S
# dat.df %>% 
#   ggplot(aes(y = stab, x = S, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lmrob",
#               formula = y~x+I(x^2)) +
#   theme_bw()
# 
# # stability ~ S
# # slight decrease in stability as size increases
# dat.df %>% filter(pca1 < -1) %>%
#   ggplot(aes(y = stab, x = S, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lmrob",
#               formula = y~x+I(x^2)) +
#   theme_bw()
# 
# # stability ~ max.tl
# dat.df %>% filter(pca1 < -1,
#                   !is.na(max.TL)) %>%
#   ggplot(aes(y = stab, x = log10(max.TL), color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
# 
# # impacted sites
# # Size
# dat.df %>% filter(pca1 > -1) %>%
#   ggplot(aes(y = stab, x = S, color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
# # max.TL
# dat.df %>% filter(pca1 > -1,
#                   !is.na(max.TL)) %>%
#   ggplot(aes(y = stab, x = log10(max.TL), color = .id)) +
#   geom_point(position = position_jitter(width = 0.3),
#              alpha = 0.1) +
#   stat_smooth(method = "lm") +
#   theme_bw()
