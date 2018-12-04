# Montoya 2009

# C = jacobian matrix

# Bi = equilibrium biomass = mi * ni
# mi = mean bodymass
# ni = numerical abundance

# aij = interaction coefficient = (mj / mi)^x

# A = matrix of aij
# D = diagonal matrix of equilibrium biomass = Dii = Bi
# C = DA
# cij = -aij * Bj
# cji (positive interactions) = eij * aij * Bi
# eij = efficiency

# food web functions from Petchey
source("FoodWebFunctions.r")
source("6a_functions.R")
source("stability_fxns_Sauve.R")
library(plyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(forcats)
library(purrr)

# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# remove portal site bc only has one taxa
Pij$Portal <- NULL
# read in data
dat <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
# sort data by avg.dw so matrices are size ordered below
dat <- dat[with(dat, order(dat$avg.dw)),]
# split into list by site
dat.list <- split(dat, list(dat$site))
dat.list$Portal <- NULL

bodymass <- map(dat.list, pull,"avg.dw")
numerical.density <- map(dat.list, pull,"density")

montoya_stability <- function(p.matr, bodymass, numerical.density,
                              s2 = 1,
                              aij_coeff = 0.66,
                              trials = 1){
# adjacency matrix
  result = list()
  for(i in 1:trials){
    A <- rm_cycle(b_trial(p.matr)) # adjacency matrix 
    #aij = 1 if col j eats row i
    A2 <- -A
    # A2 represents negative effect of col j on row i
    A2[which(A2 > t(A2))] <- -t(A2)[which(A2 > t(A2))]
    # add positive effects of row j on col i
  
    S <- length(bodymass) # number od species S
    # equilibrium biomass identity matrix
    B <- bodymass * numerical.density
    D <- matrix(0, S, S)
    diag(D) <- B
    
    # interaction coefficient matrix M
    M <- t(outer(bodymass, bodymass, FUN = "/")^aij_coeff)
    C <- D%*%M
    C <- C*A2
    pos.C.index <- which(C>0, arr.ind = TRUE)
    pos.C.strength <- C[pos.C.index]
    neg.C.index <- cbind(pos.C.index[,2], pos.C.index[,1])
    neg.C.strength <- -pos.C.strength
    C[neg.C.index] <- neg.C.strength
    C[C>0] <- C[C>0]* 0.7 # runif(n= 1, min = 0.6, max = 0.8)
    stab = stability(C, s2 = s2)
    fw_meas <- Get.web.stats(A)
    # diag.strength <- colSums(A)
    # diag.strength[diag.strength > 0] <- 1e-6
    # diag.strength[diag.strength == 0] <- -1
    # diag(C) <- diag.strength
    # e.vals <- Re(eigen(C)$values)
    # all.negative.eigen <- all(Re(eigen(C)$values) < 0)
    result[[i]] <- c(stab = stab, fw_meas)  
    }
  return(result)
}

#montoya_stability(p.matr = Pij[[1]], bodymass = dat.list[[1]]$avg.dw,
#                   numerical.density = dat.list[[1]]$density)
set.seed(3049)
stab <- pmap(list(p.matr = Pij,
                  bodymass = bodymass,
                  numerical.density = numerical.density),
             montoya_stability,
             trials = 250, s2 = 2, aij_coeff = 0.75)
# all.negative <- NULL
# for(i in 1:length(stab)){
#   for (j in 1:length(stab[[i]])){
#     out <- data.frame(site = names(stab)[i],
#                      all.negative = stab[[i]][[j]][[2]])
#     all.negative <- rbind(all.negative, out)
#   }
# }
# all.negative %>% group_by(site) %>%
#   summarize(n = n(), all.negative = sum(all.negative)) %>%
#   mutate(prop.stable = all.negative / n) %>%
#   left_join(gradient, by = "site") %>%
#   ggplot(aes(y = prop.stable, x = pca1)) +
#   geom_point()

stab <- map(stab, ldply)
stab <- ldply(stab)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
stab <- left_join(stab, gradient, by = c(".id" = "site"))

joy.stability(stab, title = "Montoya et. al 2009", scale = 3, rel_min_height = 0.007, bandwidth = NULL, xmin = -.0001, xmax = 0.003)

stab %>% filter(stab < 0.02) %>%
  plot(log10(stab + 1)~pca1, data = .)

stab %>% 
  filter(stab < 0.02) %>%
  ggplot(aes(y = stab, x = pca1))+
  geom_point()

fwpointrange(stab, y = "C", size = 5) +
  theme(axis.title = element_text(size = 30)) +
  #labs(title = "Montoya Connectance") +
  NULL
fwpointrange(stab, y = "L") +
  labs(title = "Montoya Connectance") +
  NULL

stab %>% 
  filter(stab < 0.02) %>%
  fwpointrange(y = "stab") +
  labs(title = "Montoya Connectance") +
  stat_smooth(method = "lm")

ggplot(stab, aes(y = C, x = pca1))+
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) + sd(x),
               geom = "pointrange", fatten = 7) +
  theme_bw() +
  labs(x = "Mining gradient")+
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 20)) 

ggplot(stab, aes(y = L, x = pca1))+
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x),
               fun.ymax = function(x) mean(x) + sd(x),
               geom = "pointrange", fatten = 7) +
  theme_bw() +
  labs(x = "Mining gradient")+
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 20)) 

stab %>% 
  filter(stab < 0.02) %>%
ggplot( aes(y = stab, x = pca1))+
  stat_summary(fun.y = mean,
             fun.ymin = function(x) mean(x) - (sd(x) / sqrt(length(x))),
             fun.ymax = function(x) mean(x) + (sd(x) / sqrt(length(x))),
             geom = "pointrange", fatten = 5) +
  geom_hline(aes(yintercept = 0),
             linetype = 2) +
  theme_bw() +
  labs(x = "Mining gradient",
       y = "Stability") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 20))

stab %>% 
  filter(stab < 0.02) %>%
  ggplot(aes(y = stab, x = pca1))+
  geom_point() +
  geom_smooth(method = 'gam')

ggplot(stab, aes(x = stab,
                 y = fct_reorder(.id, S),
                 height = ..density.., 
                 fill = pca1)) +
  stat_density_ridges(scale = 3, rel_min_height = 0.007,
                      bandwidth = NULL, alpha = 0.8) +
  theme_ridges(center_axis_labels = TRUE) +
  labs(y = "S", x = "stability") +
  scale_fill_distiller(palette = "Spectral") +
  theme(legend.position = "NULL")+
  coord_cartesian(xlim = c(-0.0001, 0.003
                           )) +
    labs(title = "Montoya") 

ggplot(stab, aes(y = S, x = pca1))+
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "Mining gradient",
       y = "Richness") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 20))

# save S and C parameters for Random matrices
S_C_parameters <- stab %>%
  group_by(.id) %>%
  summarise(S.mean = mean(S),
            S.sd = sd(S),
            C.mean = mean(C),
            C.sd = sd(C))

S_C_parameters <- S_C_parameters %>%
  mutate(C = C.mean,
         C.hi = C.mean + C.sd,
         C.lo = C.mean - C.sd)

S_C_parameters <- S_C_parameters %>%
  select(S.mean, C, C.hi, C.lo) %>%
  gather("C_param", "C_value", 2:4) %>%
  arrange(S.mean)

saveRDS(S_C_parameters, "data/S_C_parameters_for_random_matrices.RDS")
