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
    #A <- b_trial(matr)
    A <- matr
    A <- rm_low_tri(adjacency = A)
    xistar <- get_xistar(N = dat$no.m2,
                         M = dat$dw)
    aij <- get_aij(xistar)
    mij <- get_mij(A = A, aij = aij, xi = xistar)
    #diag(mij) <- -1
    eig1 <- Re(eigen(mij)$values)[1]
    eigmax <- max(Re(eigen(mij)$values))
    return.time = -1 / eig
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(eig1 = eig1,
                     eigmax= eigmax,
                     return.time = return.time,
                     fw_meas) 
  }
  return(result)
}


# data ####
# invertebrate biomass and abundance
invert <- readRDS("C:\\Users\\Justin\\Google Drive\\Data\\Predicting NZ Food Webs\\estimated invert bodymass.RDS")
# fish biomass and abundance
# the RDS file is a list of 4 different fish sizes
# [[2]] is the mean minimum 
fish <- readRDS("C:\\Users\\Justin\\Google Drive\\Data\\Predicting NZ Food Webs\\estimated fish bodymass and abundance.RDS")[[2]]
# rbind invert to each fish, filter out na(dw), order by dw and then split into list of sites
dw <- bind_rows(fish, invert)
dw <- dw[!is.na(dw$dw),]
dw <- dw[order(dw$dw),]
dw <- split(dw, list(dw$site))

obs.A <- readRDS("C:\\Users\\Justin\\Google Drive\\Data\\Predicting NZ Food Webs\\observed matrices matched to inferred.RDS")

match_matr <- function (A, taxa){
  # colnames(inf) have already been size-sorted
  index = intersect(taxa$taxa, rownames(A))
  A = A[index, index]
  #inf = inf[index, index]
  #list(observed = obs, inferred = inf)
  A
}

obs.A <- map2(obs.A, dw, match_matr)

results <- map2(test, dw, get_measures, trials = 1)
results <- map(results, ldply)
results <- ldply(results)



A <- obs.A$Dempsters

M <- rm_low_tri(adjacency = A)
xistar <- get_xistar(N = dw$Dempsters$no.m2,
                     M = dw$Dempsters$dw)
aij <- get_aij(xistar)
mij <- get_mij(A = M, aij = aij, xi = xistar)
eye.approximate.ReL1(mij)
