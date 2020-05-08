# 5_final_probability_matrices

# estimate link probability by multiplying niche probability (P.niche) and neutral probability (N)
# prune niche forbidden links e.g. Pomeranz et al. 2019 Methods Eco Evo

# This is based on Poisot et al 2014: 
# Pij [proportional to] P.niche(i,j) * N(i,j)
# e.g. probability of interaction is proportional to the product of niche and neutral effects
# P.niche = probability calculated by predict.niche.prob() in script #2_link_probability.R
# N = ni * nj where n is a vector of relative abundances calculated in script 4_relative_abundance_matrices.R 

library(plyr)
library(tidyverse)
library(gplots)

# functions written for this MS
source("functions/MS_functions.R")

# data ####
# read in link probabilities (P.niche)
P.niche <- readRDS("data/AMD_niche_link_probability_matrices.RDS")
# read in relative abundance matrices
N <- readRDS("data/AMD_relative_abundance_matrices.RDS")

# niche forbidden ####
# "forbidden" taxa based on morphology
# e.g. "brushing/scraping" mouthparts, 
taxa.forbid <- c("Acari", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Copepoda", "Deleatidium", "Elmidae", "Elmidae Adult", "Helicopsyche", "Hydraenidae", "Nesameletus","Oligochaetae", "Olinga","Ostracoda", "Oxyethira", "Paraleptamphopus", "Potamopyrgus", "Rakiura", "Scirtidae", "Spaniocerca", "Spaniocercoides", "Zelandobius", "Zelolessica", "Zephlebia")

# write taxa forbid as table for MS
data.frame(Taxa = taxa.forbid) %>%
  write.csv("figures/Niche_forbid_taxa.csv",
            row.names = FALSE)

# set forbidden taxa columns = 0
P.niche <- map(P.niche,
         rm_niche, # in MS_functions.R
         taxa = taxa.forbid)

# P.niche * relative abundance matrices
Pij <- map2(P.niche, N, ~.x*.y)

# rescale probabilities between 0 and 1
Pij <- map(Pij,
           scalexy, # function in MS_functions.R
           min = 0.01,
           max = 0.99)

# save final probabilities which account for
# niche, neutral, and forbidden links
saveRDS(Pij, "data/AMD_final_probability_matrices.RDS")
