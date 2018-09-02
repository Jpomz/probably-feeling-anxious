# scale interaction strengths based on row and column ranks

# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(forcats)
library(gplots)

# food web functions from Petchey
source("FoodWebFunctions.r")
# functions
source("6a_functions.R")
source("stability_fxns_Sauve.R")

# function to get fw_measures and dominant eigenvalue
get_measures <- function(matr, s2,
                         trials = 1,
                         scale.Jij = FALSE,
                         correlate.Jij = FALSE){
  result <- list()
  for(i in 1:trials){
    A <- rm_cycle(b_trial(matr))
    J <- jacobian_binary(A)
    if(scale.Jij == TRUE){
      J = scale.Jij(J)
    }
    if(correlate.Jij == TRUE){
      J = correlate.Jij(J)
    }
    stab = stability(J, s2 = s2)
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(stab = stab, fw_meas) 
  }
  return(result)
}

joy.stability <- function(data, x = "stab", xmin = -.09, xmax = 0.5,
                          title = NULL, scale = NULL,
                          bandwidth = NULL,
                          rel_min_height = NULL){
  ggplot(data, aes(x = data[[x]],
                   y = fct_reorder(.id, pca1),
                   height = ..density.., 
                   fill = pca1)) +
    stat_density_ridges(scale = scale, rel_min_height = rel_min_height,
                        bandwidth = bandwidth, alpha = 0.8) +
    theme_ridges(center_axis_labels = TRUE) +
    labs(y = "mining gradient", x = "stability") +
    scale_fill_distiller(palette = "Spectral") +
    theme(legend.position = "NULL")+
    coord_cartesian(xlim = c(xmin, xmax)) +
    if (is.null(title))
      labs(title = NULL) 
  else
    labs(title = title)
}

# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# remove portal site bc only has one taxa
Pij$Portal <- NULL

# random interaction strengths
# for reproducibility
set.seed(3049)
random.J <- map(Pij, get_measures, s2 = 1, trials = 250)
random.J <- map(random.J, ldply)
random.J <- ldply(random.J)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
random.J <- left_join(random.J, gradient, by = c(".id" = "site"))
joy.stability(random.J, title = "Random Jij",
              scale = 3, rel_min_height = 0.07,
              bandwidth = 0.007, xmin = -0.02, xmax = 0.13)
fwpointrange(random.J, y = "C") +
  labs(title = "random.J Connectance")
ggplot(random.J, aes(y = stab, x = pca1))+
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) +sd(x),
               geom = "pointrange", fatten = 5) +
  theme_bw() +
  labs(x = "Mining gradient",
       y = "Stability") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 20)) 


# scaled interaction strengths
# for reproducibility
set.seed(3049)
scaled.J <- map(Pij, get_measures, s2 = 2, trials = 250, scale.Jij = TRUE)
scaled.J <- map(scaled.J, ldply)
scaled.J <- ldply(scaled.J)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
scaled.J <- left_join(scaled.J, gradient, by = c(".id" = "site"))
joy.stability(scaled.J, title = "Body size scaled Jij", scale = 3, rel_min_height = 0.07, bandwidth = 0.007, xmin = -0.02, xmax = 0.13)
#fwpointrange(scaled.J, y = "C") +
#  labs(title = "scaled.J Connectance")


# correlated interaction strengths
# for reproducibility
set.seed(3049)
corr.J <- map(Pij, get_measures, s2 = 2, trials = 250, correlate.Jij = TRUE)
corr.J <- map(corr.J, ldply)
corr.J <- ldply(corr.J)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
corr.J <- left_join(corr.J, gradient, by = c(".id" = "site"))
joy.stability(corr.J, title = "Correlated Jij", scale = 3, rel_min_height = 0.07, bandwidth = 0.007, xmin = -0.02, xmax = 0.13)

# scaled and correlated interaction strengths
# for reproducibility
set.seed(3049)
s.c.J <- map(Pij, get_measures,
              s2 = 2, trials = 250, correlate.Jij = TRUE, scale.Jij = TRUE)
s.c.J <- map(s.c.J, ldply)
s.c.J <- ldply(s.c.J)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
s.c.J <- left_join(s.c.J, gradient, by = c(".id" = "site"))
joy.stability(s.c.J, title = "Scaled and correlated Jij",
              scale = 3, rel_min_height = 0.07,
              bandwidth = 0.007, xmin = -.02, xmax = 0.13)

# save S and C parameters for Random matrices
S_C_parameters <- random.J %>%
  group_by(.id) %>%
  summarise(S.mean = mean(S),
            S.sd = sd(S),
            C.mean = mean(C),
            C.sd = sd(C))

S_C_parameters <- S_C_parameters %>%
  mutate(C = C.mean,
         C.hi = C.mean + C.sd,
         C.lo = C.mean - C.sd)
# write table for MS
S_C_parameters %>%
  select(S.mean, C, C.hi, C.lo) %>%
  arrange(S.mean) %>%
  write.csv("figures/S_C_Values.csv",
            row.names = FALSE)

S_C_parameters <- S_C_parameters %>%
  select(S.mean, C, C.hi, C.lo) %>%
  gather("C_param", "C_value", 2:4) %>%
  arrange(S.mean)

saveRDS(S_C_parameters, "data/S_C_parameters_for_random_matrices.RDS")


stability.results <- list(random = random.J[,-1], scaled = scaled.J[,-1], corr = corr.J[,-1], s.c = s.c.J[,-1])
stability.results <- ldply(stability.results)

saveRDS(stability.results, "data/stability_results.RDS")






get_Jij <- function(matr, trials = 1, scale.Jij = FALSE, 
            correlate.Jij = FALSE, plot = FALSE, 
            density = FALSE, ...){
  result <- list()
  for(i in 1:trials){
    A <- b_trial(matr)
    A <- rm_cycle(A)
    low.tri <- sum(A[lower.tri(A)])
    prop.low.tri <- low.tri / sum(A)
    J <- jacobian_binary(A)
    result[[i]] <- data.frame(sum = low.tri,
                              percent = prop.low.tri)
    if(scale.Jij == TRUE){
      J = scale.Jij(J)
    }
    if(correlate.Jij == TRUE){
      J = correlate.Jij(J)
    }
    if(density == TRUE){
      interactions = J[J!=0]
      plot(density(interactions))
    }
    if(plot == TRUE){
    heatmap.2(J, 
              Rowv = NA, Colv = NA, scale = "none",
              trace = "none",
              dendrogram = "none", key = F, labRow = NA,
              labCol = NA)
    }
  }
  return(result)
}
# Smanually scale colors
# col = colorpanel(n = 50, low = "blue", mid = "beige", high = "red")
A <- rm_cycle(b_trial(Pij$Lankey))
heatmap.2(A, 
          Rowv = NA, Colv = NA, scale = "none",
          trace = "none",
          dendrogram = "none", key = F, labRow = NA,
          labCol = NA)
J <- jacobian_binary(A)
heatmap.2(correlate.Jij(scale.Jij(J)), 
          Rowv = NA, Colv = NA, scale = "none",
          trace = "none",
          dendrogram = "none", key = F, labRow = NA,
          labCol = NA)

set.seed(3049)
low <- map(Pij, get_Jij, trials = 250, rm_cycle = TRUE)
low <- map(low, ldply)
low <- ldply(low)

ggplot(low, aes(y = percent, x = .id))+
  geom_point() +
  coord_cartesian(ylim = c(0, 0.6))

# heatmap.2(Pij[[1]], 
#   Rowv = NA, Colv = NA, scale = "none", trace = "none",
#   dendrogram = "none", key = F, labRow = NA,
#   labCol = NA)




lm.list <- map(stability.results, ~lm(stab~pca1, data = .x))

lm.list %>%
  map(~summary(.x))

stability.results %>%
  map(~lm(log10(stab + 1.1)~pca1, data = .x)) %>%
  map(~summary(.x))


x <- seq(1e-6, 1e-5, length.out = 100)
expand.grid(x,x) %>% 
  filter( Var1 > Var2) %>% 
  mutate(aij = (Var1 / Var2)**0.75) %>% 
  plot(aij ~ Var1, data = .)


# plots for powerpoint outline
heatmap.2(J, 
          Rowv = NA, Colv = NA, scale = "none",
          trace = "none",
          dendrogram = "none", key = F, labRow = NA,
          labCol = NA, col = colorpanel(n = 50, low ="turquoise1" , mid = "lavenderblush", high = "red"))
