# interaction strengths
# sensu Montoya et al. 2005

# interaction strength ~
# aij = log10((Mc / Mr)^0.75)

library(plyr)
library(tidyverse)
# bernoulli trials
dat <- readRDS("data/500 bernouli trials.rds")
amd <- readRDS("data/AMD fish invert dw abundance.RDS")
amd.list <- split(amd, list(amd$site))
amd.list <- map(amd.list, arrange, avg.dw)
dw.list <- map(amd.list, pull, "avg.dw")

get_intxn_strength <- function(web, dw, n){
  pairs <- which(web > 0, arr.ind = TRUE)
  consumer <- dw[pairs[, 2]]
  resource <- dw[pairs[, 1]]
  intxn <- (consumer / resource)**0.75
  return(intxn)
}


intxns <- NULL
for (site in 1:length(dat)){
  #c(1:11, 13:15, 17:20, 22:25)
  out <- NULL
  for (i in 1:length(dat[[site]])){
    if (sum(dat[[site]][[i]]) <= 1){
      next
      } else {
      out[[i]] <- get_intxn_strength(
        web = dat[[site]][[i]],
        dw = dw.list[[site]])
    } 
  }
  intxns[[site]] <- out
}

out <- NULL
for (i in 1:length(dat[[12]])){
  if (sum(dat[[12]][[i]]) <= 1){#if (sum(dat[[site]][[i]]) == 0)
    next #out[[i]] <- NA
    } else {
    out[[i]] <- get_intxn_strength(
      web = dat[[12]][[i]],
      dw = dw.list[[12]])
  } 
}



get_intxn_strength(
  web = dat[[12]][[3]],
  dw = dw.list[[12]])
