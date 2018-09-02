# 6_food_web_calculations
# jfpomeranz@gmail.com
# June 2018

# calculate dominant eigenvalue and food web measures
# from probability matrices and data on taxa bodymass and abundance

# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)

# food web functions from Petchey
source("FoodWebFunctions.r")
# functions
source("6a_functions.R")

# function to get fw_measures and dominant eigenvalue
get_measures <- function(matr, dat, trials = 1){
  result <- list()
  for(i in 1:trials){
    A <- b_trial(matr)
    A <- rm_low_tri(adjacency = A)
    xistar <- get_xistar(N = dat$density, M = dat$avg.dw)
    aij <- get_aij(xistar)
    mij <- get_mij(A = A, aij = aij, xi = xistar)
    #diag(mij) <- -1
    eig <- Re(eigen(mij)$values)[1]
    #eig <- max(Re(eigen(mij)$values))
    return.time = -1 / eig
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(dom.eig = eig,
                     #return.time = return.time,
                     fw_meas) 
  }
  return(result)
}

# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# taxa data
taxa.dat <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
# order taxa.dat by size to match Pij
dat <- taxa.dat[order(taxa.dat$avg.dw),]
# split taxa data into list by site
dat <- split(taxa.dat, list(taxa.dat$site))

# for reproducibility
set.seed(3049)

# remove portal site bc only has one taxa
Pij[[21]] <- NULL
dat[[21]] <- NULL

results <- map2(Pij, dat, get_measures, trials = 50)
results <- map(results, ldply)
results <- ldply(results)
gradient <- readRDS("data/pca_axis_26_sites.RDS")
results <- left_join(results, gradient, by = c(".id" = "site"))

saveRDS(results, "data/AMD_fw_measures.RDS")
