# Stability statistics
library(dplyr)
library(purrr)
library(ggplot2)
library(robustbase)
library(lmtest)
library(forcats)

dat <- readRDS("data/random_network_stability_results.RDS") %>%
  split(list(.$.id))
# heteroskedastic?
dat %>% map(~bptest(stab~S, data = .x))       
# all are

# linear or quadratic?
dat %>% map(~anova(lmrob(stab~S, data = .x,
                         control = lmrob.control(
                           maxit.scale = 700,
                           rel.tol = 1e-6)),
                   lmrob(stab~S +I(S^2), data = .x,
                         control = lmrob.control(
                           maxit.scale = 700,
                           rel.tol = 1e-6)),
                   test = "Wald"))
# only s.c is quadratic

dat %>% map(~summary(lmrob(stab~S +I(S^2), data = .x,
                           control = lmrob.control(
                             maxit.scale = 700,
                             rel.tol = 1e-6))))

dat$random %>%
  lmrob(stab~S, data = .) %>%
  summary()

fig_fun <- function(data, mainlab,
                    x = "S", y = "stab"){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_discrete(h.start = 80) +
    labs(x = "Network size",
         y = "stability metric",
         title = mainlab)+
    stat_smooth(#method = "lm",
       method = "lmrob",
                 #formula = y~x+I(x^2),
                color = "black") 
}
fig_fun(data = dat$random, mainlab = "Random structure and interaction strength")
ggsave("figures/s~R_str_R_strenth.png")

dat %>%
  map(~fig_fun(data = .x, mainlab = paste(.x[1,1])))
