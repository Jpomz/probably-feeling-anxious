# Different interaction strength estimates
# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(forcats)
library(gplots)

# food web functions from Petchey
source("functions/FoodWebFunctions.r")
# functions written for this MS
source("functions/MS_functions.R")
# stability functions from Sauve 2016
# if using, please cite original publication:
# Sauve, A. M. C., Thébault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1
source("functions/stability_fxns_Sauve.R")


# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# gradient
gradient <- readRDS("data/pca_axis_26_sites.RDS")


# random interaction strengths ####
# set seed for reproducibility
set.seed(3049)
random.J <- ldply(
  map(
    map(
      Pij,
      get_measures, 
      s2 = 2,
      trials = 250),
    ldply)
)
# join gradient + random.J
random.J <- left_join(random.J, gradient, by = c(".id" = "site"))

# plots
# ggplot(random.J, aes(x = stab,
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.05,
#                        scale = 8,
#                        alpha = 0.45) +
#   scale_fill_manual(values = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme_bw() +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15)) +
#   labs(title = "Random Jij",
#        y = "Site",
#        x = "s") +
#   coord_cartesian(xlim = c(-.009, .11)) +
#   NULL

# ggplot(random.J, aes(x = C,
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.05,
#                       scale = 1,
#                       alpha = 0.45) +
#   scale_fill_manual(values = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme_bw() +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15)) +
#   labs(y = "Site",
#        x = "C") +
#   #coord_cartesian(xlim = c(-.009, .11)) +
#   NULL
# 
# ggplot(random.J, aes(x = log10(L),
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.01,
#                       scale = 10,
#                       alpha = 0.45) +
#   scale_fill_manual(values = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme_bw() +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15)) +
#   labs(y = "Site",
#        x = "L") +
#   #coord_cartesian(xlim = c(-.009, .11)) +
#   NULL




# scaled interaction strengths ####
# for reproducibility
set.seed(3049)
scaled.J <- ldply(
  map(
    map(
      Pij,
      get_measures,
      s2 = 2,
      trials = 250,
      scale.Jij = TRUE),
    ldply)
  )

scaled.J <- left_join(scaled.J, gradient, by = c(".id" = "site"))

# plots
# ggplot(scaled.J, aes(x = stab,
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.05,
#                        scale = 8,
#                        alpha = 0.45) +
#   theme_bw() +  
#   scale_fill_manual(values = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15)) +
#   labs(title = "Scaled Jij",
#        y = "Site",
#        x = "s") +
#   coord_cartesian(xlim = c(-.004, .05)) +
#   NULL

# correlated interaction strengths ####
# for reproducibility
set.seed(3049)
corr.J <- ldply(
  map(
    map(
      Pij,
      get_measures, 
      s2 = 2, 
      trials = 250, 
      correlate.Jij = TRUE),
    ldply)
)
corr.J <- left_join(corr.J, gradient, by = c(".id" = "site"))
#plots
# ggplot(corr.J, aes(x = stab,
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.05,
#                       scale = 8,
#                       alpha = 0.45) +
#   theme_bw() +  
#   scale_fill_manual(values = colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 15)) +
#   labs(title = "Scaled Jij",
#        y = "Site",
#        x = "s") +
#   #coord_cartesian(xlim = c(-.004, .05)) +
#   NULL

# scaled and correlated interaction strengths ####
# for reproducibility
set.seed(3049)
s.c.J <- ldply(
  map(
    map(
      Pij,
      get_measures,
      s2 = 2,
      trials = 250, 
      correlate.Jij = TRUE,
      scale.Jij = TRUE),
    ldply),
)

s.c.J <- left_join(s.c.J, gradient, by = c(".id" = "site"))

# ggplot(s.c.J, aes(x = stab,
#                      y = fct_reorder(.id, pca1),
#                      fill = fct_reorder(.id, pca1))) +
#   geom_density_ridges(rel_min_height = 0.05,
#                        scale = 3,
#                        alpha = 0.5) +
#   theme_bw() +
#   scale_fill_manual(
#     values = colorRampPalette(rev(
#       brewer.pal(n = 11, name ="Spectral")))(25)) +
#   theme(legend.position = "NULL",
#         axis.title = element_text(size = 10),
#         axis.text = element_text(size = 10)) +
#   labs(title = "Scaled & Correlated Jij",
#        y = "Site",
#        x = "s") +
#   coord_cartesian(xlim = c(-.009, .11)) +
#   NULL


# save S and C parameters for Random matrices
# S_C_parameters <- random.J %>%
#   group_by(.id) %>%
#   summarise(S.mean = mean(S),
#             S.sd = sd(S),
#             C.mean = mean(C),
#             C.sd = sd(C))
# 
# S_C_parameters <- S_C_parameters %>%
#   mutate(C = C.mean,
#          C.hi = C.mean + C.sd,
#          C.lo = C.mean - C.sd)
# write table for MS
# S_C_parameters %>%
#   select(S.mean, C, C.hi, C.lo) %>%
#   arrange(S.mean) %>%
#   write.csv("figures/S_C_Values.csv",
#             row.names = FALSE)
# 
# S_C_parameters <- S_C_parameters %>%
#   select(S.mean, C, C.hi, C.lo) %>%
#   gather("C_param", "C_value", 2:4) %>%
#   arrange(S.mean)
# 
# saveRDS(S_C_parameters, "data/S_C_parameters_for_random_matrices.RDS")

# combine all stability results into a list. 
# removing column 1 (".id") from each dataframe. This information is redundant with "Site" column
stability.results <- list(random = random.J[,-1],
                          scaled = scaled.J[,-1],
                          corr = corr.J[,-1],
                          s.c = s.c.J[,-1])
stability.results <- ldply(stability.results)

saveRDS(stability.results, "data/stability_results.RDS")



# below this is trash? ####
# 
# 
# get_Jij <- function(matr, trials = 1, scale.Jij = FALSE, 
#             correlate.Jij = FALSE, plot = FALSE, 
#             density = FALSE, ...){
#   result <- list()
#   for(i in 1:trials){
#     A <- b_trial(matr)
#     A <- rm_cycle(A)
#     low.tri <- sum(A[lower.tri(A)])
#     prop.low.tri <- low.tri / sum(A)
#     J <- jacobian_binary(A)
#     result[[i]] <- data.frame(sum = low.tri,
#                               percent = prop.low.tri)
#     if(scale.Jij == TRUE){
#       J = scale.Jij(J)
#     }
#     if(correlate.Jij == TRUE){
#       J = correlate.Jij(J)
#     }
#     if(density == TRUE){
#       interactions = J[J!=0]
#       plot(density(interactions))
#     }
#     if(plot == TRUE){
#     heatmap.2(J, 
#               Rowv = NA, Colv = NA, scale = "none",
#               trace = "none",
#               dendrogram = "none", key = F, labRow = NA,
#               labCol = NA)
#     }
#   }
#   return(result)
# }
# # manually scale colors
# # col = colorpanel(n = 50, low = "blue", mid = "beige", high = "red")
# A <- rm_cycle(b_trial(Pij$Lankey))
# heatmap.2(A, 
#           Rowv = NA, Colv = NA, scale = "none",
#           trace = "none",
#           dendrogram = "none", key = F, labRow = NA,
#           labCol = NA)
# J <- jacobian_binary(A)
# heatmap.2(correlate.Jij(scale.Jij(J)), 
#           Rowv = NA, Colv = NA, scale = "none",
#           trace = "none",
#           dendrogram = "none", key = F, labRow = NA,
#           labCol = NA)
# 
# set.seed(3049)
# low <- map(Pij, get_Jij, trials = 250, rm_cycle = TRUE)
# low <- map(low, ldply)
# low <- ldply(low)
# 
# ggplot(low, aes(y = percent, x = .id))+
#   geom_point() +
#   coord_cartesian(ylim = c(0, 0.6))
# 
# # heatmap.2(Pij[[1]], 
# #   Rowv = NA, Colv = NA, scale = "none", trace = "none",
# #   dendrogram = "none", key = F, labRow = NA,
# #   labCol = NA)
# 
# 
# 
# 
# lm.list <- map(stability.results, ~lm(stab~pca1, data = .x))
# 
# lm.list %>%
#   map(~summary(.x))
# 
# stability.results %>%
#   map(~lm(log10(stab + 1.1)~pca1, data = .x)) %>%
#   map(~summary(.x))
# 
# 
# x <- seq(1e-6, 1e-5, length.out = 100)
# expand.grid(x,x) %>% 
#   filter( Var1 > Var2) %>% 
#   mutate(aij = (Var1 / Var2)**0.75) %>% 
#   plot(aij ~ Var1, data = .)
# 
# 
# # plots for powerpoint outline
# heatmap.2(J, 
#           Rowv = NA, Colv = NA, scale = "none",
#           trace = "none",
#           dendrogram = "none", key = F, labRow = NA,
#           labCol = NA, col = colorpanel(n = 50, low ="turquoise1" , mid = "lavenderblush", high = "red"))
