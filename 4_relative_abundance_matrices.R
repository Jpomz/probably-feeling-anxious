# 3_relative_abundance_matrices
# june 2018
# jfpomeranz@gmai.com

# 1) calculate relative abundance matrices for AMD sites
# 2) "correct" fish relative abundances as in 
# Pomeranz et al. in review Methods in Eco and Evo
# 3) rescale relative abundances to be between 0.5 and 1

library(plyr)
library(purrr)

# functions written for this MS
source("functions/MS_functions.R")

# read in data
dat <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
# sort data by avg.dw so matrices are size ordered below
dat <- dat[with(dat, order(dat$avg.dw)),]
# split into list by site
dat.list <- split(dat, list(dat$site))

# 1) relative abundance ####
# calculate relative abundance matrices
rel.ab <- llply(dat.list, function (x){
  get_rel_ab(vec = x$rel.ab, taxa = x$taxa)
})

# 2) correct fish ####
# vector of fish taxa
taxa.fish <- c("Salmo.trutta", "Anguilla.australis",
               "Galaxias.fasciatus", "Galaxias.maculatus",
               "Gobiomorphus.huttoni")
rel.ab <- map(rel.ab, correct_fish_rel_ab, fish = taxa.fish)

# 3) rescale ####
rel.ab <- map(rel.ab, scalexy, min = 0.5, max = 1)

saveRDS(rel.ab, "data/AMD_relative_abundance_matrices.RDS")

