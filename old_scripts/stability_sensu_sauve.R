# 6_food_web_calculations
# jfpomeranz@gmail.com
# July 2018

# Estimate stability sensu and functions of Sauve et al. 2016
# using random interaction strengths drawn from rnorm(mu = 1, sd = 0.1)


# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(forcats)

# food web functions from Petchey
source("FoodWebFunctions.r")
# functions
source("6a_functions.R")
source("stability_fxns_Sauve.R")

# function to get fw_measures and dominant eigenvalue
get_measures <- function(matr, s2, trials = 1){
  result <- list()
  for(i in 1:trials){
    A <- b_trial(matr)
    A <- rm_low_tri(adjacency = A)
    J <- jacobian_binary(A)
    stab = stability(J, s2 = s2)
    #fw_meas <- Get.web.stats(A)
    result[[i]] <- c(stab = stab)#, fw_meas) 
  }
  return(result)
}

# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# remove portal site bc only has one taxa
Pij[[21]] <- NULL

# taxa data
# taxa.dat <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
# # order taxa.dat by size to match Pij
# dat <- taxa.dat[order(taxa.dat$avg.dw),]
# # split taxa data into list by site
# dat <- split(taxa.dat, list(taxa.dat$site))
#dat[[21]] <- NULL

# for reproducibility
set.seed(3049)
results <- map(Pij, get_measures, s2 = 1, trials = 250)
results <- map(results, ldply)
results <- ldply(results)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
results <- left_join(results, gradient, by = c(".id" = "site"))
plot(stab~pca1, dat = results)
ggplot(results, aes(y = stab, x = pca1))+
  geom_point() +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) + sd(x),
               geom = "pointrange",
               color = "red") +
  stat_summary(fun.y = median,
               geom = "point",
               color = "green") +
  geom_smooth(method = "lm")

ggplot(results, aes(y = fct_reorder(.id, pca1),
                    x = stab,
                    #group = pca1,
                    height = ..density..,
                    fill = pca1)) +
  stat_density_ridges(scale = 3,
                      rel_min_height = 0.05,
                      alpha = 0.8) +
                      #quantile_lines = TRUE,
                      #quantiles = 2) + 
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "mining gradient", 
       x = "stability") +
  scale_fill_distiller(palette = "Spectral") +
  theme(legend.position = "NULL")

# ggplot(results, aes(y = fct_reorder(.id, pca1),
#                     x = stab,
#                     #group = pca1,
#                     height = ..density..,
#                     fill = pca1)) +
#   geom_density_ridges(scale = 5,
#                       #rel_min_height = 0.01,
#                       alpha = 0.5,
#                       stat = "density") + 
#   theme_ridges(center_axis_labels = TRUE) +
#   labs(y = "mining gradient", 
#        x = "stability") +
#   scale_fill_distiller(palette = "Spectral") +
#   theme(legend.position = "NULL")
# 
# results %>%
#   group_by(.id) %>%
#   do(ggplot2:::compute_density(.$stab, NULL)) %>%
#   rename(stab = x) -> stab_densities
# head(stab_densities)
# ggplot(stab_densities, aes(x = stab, y = .id, height = density)) +
#   geom_density_ridges(stat = "identity")

saveRDS(results, "data/AMD_fw_measures.RDS")