# random structure

library(tidyr)
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggridges)

source("FoodWebFunctions.r")
source("6a_functions.R")
source("stability_fxns_Sauve.R")

S_C_parameters <- readRDS("data/S_C_parameters_for_random_matrices.RDS")
S <- S_C_parameters$S.mean
C <- S_C_parameters$C_value


random_A <- function(S, C, ...){
  A <- 0
  while(sum(A) == 0){
    A <- matrix(rbinom(n = S^2, prob = C, size = 1), S, S)
    diag(A) <- 0
    rm_cycle(A)
  }
  return(A)
}

# test_function <- function(x, y){
#   A <- (x + y) / y
#   A
# }


# function to get fw_measures and dominant eigenvalue
get_random_measures <- function(S,
                                C,
                                s2,
                                trials = 1,
                                scale.Jij = FALSE,
                                correlate.Jij = FALSE){
  result <- list()
  for(i in 1:trials){
    A <- random_A(S, C)
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

# test <- get_random_measures(S = S, C = C, trials = 1, s2 = 1)
# random structure and strength ####
set.seed(3049)
r_structure_r_strength <- map2(S, C,
                               get_random_measures,
                               s2 = 2,
                               trials = 100)
names(r_structure_r_strength) <- S
r_structure_r_strength <- map(r_structure_r_strength, ldply)
r_structure_r_strength <- ldply(r_structure_r_strength)

# random structure scaled strength ####
set.seed(3049)
r_structure_scaled <- map2(S, C,
                           get_random_measures,
                           s2 = 2,
                           trials = 100,
                           scale.Jij = TRUE)
names(r_structure_scaled) <- S
r_structure_scaled <- map(r_structure_scaled, ldply)
r_structure_scaled <- ldply(r_structure_scaled)

# random structure correlated strength ####
set.seed(3049)
r_structure_corr <- map2(S, C,
                           get_random_measures,
                           s2 = 2, trials = 100,
                           correlate.Jij = TRUE)
names(r_structure_corr) <- S
r_structure_corr <- map(r_structure_corr, ldply)
r_structure_corr <- ldply(r_structure_corr)

# random structure scaled and correlated strength ####
set.seed(3049)
r_structure_s_c <- map2(S, C,
                           get_random_measures,
                           s2 = 2, trials = 100,
                           scale.Jij = TRUE,
                           correlate.Jij = TRUE)
names(r_structure_s_c) <- S
r_structure_s_c <- map(r_structure_s_c, ldply)
r_structure_s_c <- ldply(r_structure_s_c)

# combine random results in a list
rand.stab.results <- list(
  random = r_structure_r_strength[,-1], 
  scaled = r_structure_scaled[,-1],
  corr = r_structure_corr[,-1],
  s.c = r_structure_s_c[,-1])
rand.stab.results <- ldply(rand.stab.results)
saveRDS(rand.stab.results, "data/random_network_stability_results.RDS")


ggplot(r_structure_r_strength, aes(y = stab, x = S)) +
  geom_point()+
  stat_smooth(method = "lm") +
  theme_bw()

ggplot(r_structure_r_strength, aes(x = stab,
                    y = as.factor(S),
                    height = ..density..)) +
  stat_density_ridges(scale = 3, rel_min_height = 0.007,
                      bandwidth = NULL, alpha = 0.8) +
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "S", x = "stability") +
  theme(legend.position = "NULL")+
  coord_cartesian(xlim = c(-0.025, 2)) +
  labs(title = "Random")
ggplot(r_structure_scaled, aes(x = stab,
                                   y = as.factor(S),
                                   height = ..density..)) +
  stat_density_ridges(scale = 3, rel_min_height = 0.007,
                      bandwidth = NULL, alpha = 0.8) +
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "S", x = "stability") +
  theme(legend.position = "NULL")+
  coord_cartesian(xlim = c(-0.025, 2)) +
  labs(title = "R_scaled")
ggplot(r_structure_corr, aes(x = stab,
                                   y = as.factor(S),
                                   height = ..density..)) +
  stat_density_ridges(scale = 3, rel_min_height = 0.007,
                      bandwidth = NULL, alpha = 0.8) +
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "S", x = "stability") +
  theme(legend.position = "NULL")+
  coord_cartesian(xlim = c(-0.025, 2)) +
  labs(title = "R_corr")

ggplot(r_structure_s_c, aes(x = stab,
                                   y = as.factor(S),
                                   height = ..density..)) +
  stat_density_ridges(scale = 3, rel_min_height = 0.007,
                      bandwidth = NULL, alpha = 0.8) +
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "S", x = "stability") +
  theme(legend.position = "NULL")+
  coord_cartesian(xlim = c(-0.025, 2)) +
  labs(title = "R_s_c")

fwpointrange(r_structure_r_strength, y = "C", x = "S") +
  labs(title = "Random structure and strength", x = "S")

r_structure_r_strength %>%
  group_by(S) %>%
  summarize(mean.C = mean(C), sd.C = sd(C)) %>%
  mutate(lo = mean.C - sd.C, hi = mean.C + sd.C)

r_structure_r_strength %>%
  group_by(S) %>%
  summarize(mean.stab = mean(stab), sd.stab = sd(stab))
stab %>%
  group_by(S) %>%
  summarize(mean.stab = mean(stab), sd.stab = sd(stab))

dnorm(0.000623, mean = 0.0405, sd = 0.0169)
plot(density(rnorm(1000, mean = 0.0405, sd = 0.0169)))

pnorm(0.000623, mean = 0.0405, sd = 0.0169)

r_structure_r_strength %>%
  filter(S == 37) %>%
  pull(stab) %>%
  density() %>%
  plot

