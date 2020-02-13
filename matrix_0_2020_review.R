# do the number of 0's affect stability?
source("functions/MS_functions.R")
source("functions/stability_fxns_Sauve.R")
source("functions/FoodWebFunctions.r")

library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)

# read in modelling data
# for graphs at bottom 
dat <- readRDS("data/stability_results.RDS")
dat <- dat %>% 
  # select one model output e.g. random
  filter(.id == "random")

# read in probability matrices
Pij <- readRDS( "data/AMD_final_probability_matrices.RDS")
# 
mat41 <- Pij$Lankey
mat27 <- Pij$Kiwi
# # mat10 <- Pij$cbdale
mat5 <- Pij$Packtr


# stab_0 <- function(P, s2 = 2, trials = 10){
#   result <- matrix(0, nrow = trials, ncol = 2)
#   for (i in 1:trials){
#     A <- rm_cycle(b_trial(P))
#     J <- jacobian_binary(A)
#     num_0 <- length(J[J==0])/ncol(J)**2
#     stab <- stability(J, s2 = s2)
#     result[i,] = c(num_0, stab)
#   }
#   result
# }
# # 
# # 
# # # A <- matrix(0, 10, 10)
# # 
# # # get_measures(mat41, s2 = 2)
# # # b_trial(mat5)
# # 
# # 
# stab_mat41 <- stab_0(mat41, trials = 1000)
# plot(stab_mat41)
# # stab_mat27 <- stab_0(mat27, trials = 1000)
# # plot(stab_mat27)
# # stab_mat10 <- stab_0(mat10, trials = 1000)
# # plot(stab_mat10)
# stab_mat5 <- stab_0(mat5, trials = 1000)
# plot(stab_mat5)
# # 
# # plot(rbind(stab_mat41, stab_mat27, stab_mat10, stab_mat5))
# # # 
# # # x <- matrix(1, 2, 2)
# # # y <- matrix(2, 2, 2)
# # # rbind(x, y)
# # 
# # 
# stab.0 <- ldply(map(
#   Pij,
#   stab_0,
#   trials = 1000))
# 
# names(stab.0) <- c("site", "num_0", "stab")
# 
# ggplot(stab.0,
#        aes(x = num_0,
#            y = stab))+
#   geom_point(#aes(color = site),
#  alpha = 0.1) +
#   stat_smooth(method = "lm", formula = y~x) +
#   theme_bw()
# # 
# summary(lm(stab ~ num_0, data = stab.0))
# 
# Pij <- readRDS( "data/AMD_final_probability_matrices.RDS")
# P <- Pij$cbdale
# A <- rm_cycle(b_trial(P))
# 
# J <- jacobian_binary(A)
# s2 = 2
# 
# stability <- function(J, s2){
#   # J is the community matrix  (i.e. the jacobian matrix, see functions jacobian_binary() or jacobian_quant()) with a null diagonal, dim(J) = S*S
#   # s2 is an arbitrary parameter such that the real part of the greatest eigenvalue of J-s2*I is negative
#   
#   test1 <- (dim(J)[1] == dim(J)[2]) # Is J a square matrix?
#   test2 <- FALSE
#   if (test1 == TRUE){
#     S <- dim(J)[1]
#     test2 <- which(diag(J) != vector("numeric", S)) # Does J have a null diagonal?
#   }
#   
#   if ((test1 == TRUE) & (length(test2) == 0)){ # if J is a square matrix with a null diagonal
#     S <- dim(J)[1] # S is the number of species in the network
#     s1 <- 0
#     I <- diag(S)
#     E1 <- max(Re(eigen(J-s1*I, only.values = T)$values))
#     E2 <- max(Re(eigen(J-s2*I, only.values = T)$values))
#     if ((E1>=0) & (E2<0)){ # if s2 is well chosen and the system is not already stable
#       while ((s2-s1)>=10^-4){
#         stab <-(s1+s2)/2
#         E1 <-max(Re(eigen(J-stab*I, only.values = T)$values))
#         if (E1>=0){
#           s1<-stab
#         }
#         else {
#           s2<-stab
#         }
#       }
#       return(stab)
#     }
#     
#     if (E1<0){
#       # stop("J corresponds to a stable system.")
#       stab <- 0
#       return(stab)
#     }
#     if (E2>=0){
#       # stop("s2 is not high enough.")
#       stab <- NA
#       return(stab)
#     }
#   }
#   else { # if J is not a square matrix with a null diagonal
#     if (test1 == FALSE){
#       #stop("J is not a square matrix.")
#       stab <- "J is not square"
#       return(stab)
#     }
#     if (length(test2) > 0){
#       #stop("J does not have a null diagonal.")
#       stab <- "diag not null"
#       return(stab)
#     }
#   }
# }
# 
# 
# 
# A_5 <- rm_cycle(b_trial(mat5))
# J_5 <- jacobian_binary(A_5)
# stability(J_5, s2 = 2)
# rowSums(abs(J_5))
# data.frame(s = 6.1e-05, rows = rowSums(J_5))
# 
# 
# dat <- readRDS("data/stability_results.RDS")
# dat <- dat %>% filter(.id == "random")
# 
# plot(stab~C, data = dat)
# plot(stab~S, data = dat %>%
#        group_by(stab, S) %>%
#        distinct())
# plot(stab~L, data = dat)
# dat %>% mutate(L_S = (L / S)) %>%
#   plot(stab~L_S, data = .)
# 

# Matrices
# hold S, vary 0's
vary_0 <- function(S){
  result <- NULL
  m <- matrix(1, S, S)
  m <- rm_cycle(m)
  counter <- 1
  while(sum(m)>0){
    num_0 <- length(m[m==0])
    prop_0 <- num_0 / nrow(m)**2
    J <- jacobian_binary(m)
    stab <- stability(J, s2 = 2)
    index <- which(m == 1, arr.ind = TRUE)
    if(nrow(index)==1){
      break
    }else{
      random <- runif(n = nrow(index))
      index <- cbind(index, random)
      # numeric?
      index <- index[order(index[,3]),]
      m[index[1,1], index[1,2]] = 0 
      result[[counter]] <- c(num_0 = num_0,
                             prop_0 = prop_0,
                             stab = stab,
                             S = S)
      counter <- counter + 1
      }
    }
  return(result)
}
S <- c(5, 10, 15, 20, 25, 30, 35, 40)
S.vary.0 <- ldply(
  flatten(
    flatten(
      replicate(n = 10, 
                map(S, vary_0),
                simplify = FALSE))))

ggplot(data = S.vary.0, aes(x = prop_0, y = stab, color = as.factor(S))) +
  geom_point(alpha = 0.2) +
  #geom_jitter(height = NULL, width = 0.5, alpha = 0.5)+ 
  stat_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Random Matrices")

x <- replicate(n = 2, vary_0(S = 40), simplify = FALSE)
y <- flatten(x)
z <- ldply(y)
plot(stab~ prop_0, data = z)


x <- replicate(n = 5, vary_0(nrow = 40, ncol = 40), simplify = FALSE)



x <- replicate(n = 50, vary_0(nrow = 5, ncol = 5), simplify = FALSE)
y <- flatten(x)
z <- ldply(y)
plot(z)

vary_0_P <- function(P){
  result <- NULL
  m <- rm_cycle(b_trial(P))
  counter <- 1
  while(sum(m)>0){
    num_0 <- length(m[m==0])
    prop_0 <- num_0 / nrow(m)**2
    J <- jacobian_binary(m)
    stab <- stability(J, s2 = 2)
    index <- which(m == 1, arr.ind = TRUE)
    if(nrow(index)==1){
      break
    }else{
      random <- runif(n = nrow(index))
      index <- cbind(index, random)
      # numeric?
      index <- index[order(index[,3]),]
      m[index[1,1], index[1,2]] = 0 
      result[[counter]] <- c(num_0 = num_0, prop_0 = prop_0, stab = stab)
      counter <- counter + 1
    }
  }
  return(result)
}

x <- replicate(n = 5, vary_0_P(P = mat41), simplify = FALSE)
y <- flatten(x)
z <- ldply(y)
plot(stab~prop_0, data = z)

x <- replicate(n = 500, vary_0_P(P = mat5), simplify = FALSE)
y <- flatten(x)
z <- ldply(y)
plot(stab~prop_0, data = z)

x <- replicate(n = 10, vary_0_P(P = mat27), simplify = FALSE)
y <- flatten(x)
z <- ldply(y)
plot(stab~prop_0, data = z)
plot(stab~num_0, data = z)


# hold proportion of 0's constant, vary S
vary_S <- function(S, C){
  sim_C <- 0
  counter <- 1
  while(sim_C <= 0.9*C | sim_C>= 1.1*C){
    m <- matrix(runif(S^2, min = 0, max = 1), nrow = S, ncol = S)
    m[which(m <= 1.1 * C)] = 1
    m[which(m != 1)] = 0
    m <- rm_cycle(m)
    sim_C <- sum(m) / S^2
    counter <- counter + 1
    if(counter== 50){
      break
    }
  }
  J <- jacobian_binary(m)
  stab <- stability(J, s2 = 2)
  return(c(sim_C = sim_C, stab = stab, S = S, C_goal = C, counter = counter))
}

S <- c(5, 10, 15, 20, 25, 30, 35, 40)
# C <- seq(0.1, 0.4, length.out = length(S))
# 
# sim.list <- data.frame(expand.grid(S = S, C = C))

C.0.2 <- ldply(
  flatten(
    replicate(n = 100, 
              map(S, vary_S, C = 0.2),
              simplify = FALSE)))

C.0.3 <- ldply(
  flatten(
    replicate(n = 100, 
              map(S, vary_S, C = 0.3),
              simplify = FALSE)))

C.0.4 <- ldply(
  flatten(
    replicate(n = 100, 
              map(S, vary_S, C = 0.4),
              simplify = FALSE)))

C.0.1 <- ldply(
  flatten(
    replicate(n = 100, 
              map(S, vary_S, C = 0.1),
              simplify = FALSE)))

S.simulation <- rbind(C.0.1, C.0.2, C.0.3, C.0.4)
dim(S.simulation)
S.simulation <- S.simulation %>%
  group_by(C_goal)%>%
  mutate(med_C = median(sim_C)) %>%
  filter(sim_C >= med_C)

ggplot(data = S.simulation, aes(x = S, y = stab, color = as.factor(C_goal))) +
#  geom_point(position = "jitter") +
  geom_jitter(height = NULL, width = 0.5, alpha = 0.5)+ 
  stat_smooth(se = FALSE) +
  theme_bw() +
  ggtitle("Random Matrices")




# empirical / modelled food webs
#Need to load "dat"
# prop_0 ~ by s?
dat <- dat %>%
  mutate(prop_0 = (S^2 - L) / S^2)

ggplot(dat, aes(x = S, y = prop_0))+
  geom_point() +
  #stat_smooth(method = "lm") +
  stat_quantile(quantiles = c(0.1, 0.5, 0.9)) +
  NULL

ggplot(dat, aes(x = S, y = C))+
  geom_point() +
  #stat_smooth(method = "lm") +
  #stat_quantile(quantiles = c(0.1, 0.5, 0.9))
  NULL

