# 5_final_probability_matrices

# estimate link probability by multiplying niche probability (P.niche) and neutral probability (N)
# also prune niche forbidden links
# from Poisot et al 2014: 
# Pij [proportional to] P.niche(i,j) * N(i,j)
# P.niche = probability calculated by predict.niche.prob() from 
# Bartomeus, in script #2_link_probability.R
# N ~ ni * nj where n is a vector of relative abundances, rescaled from
# 0.5 to 1. 

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
# specialized filtering, Austrosimulium's cephalic fan
taxa.forbid <- c("Acari", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Copepoda", "Deleatidium", "Elmidae", "Elmidae Adult", "Helicopsyche", "Hydraenidae", "Nesameletus","Oligochaetae", "Olinga","Ostracoda", "Oxyethira", "Paraleptamphopus", "Potamopyrgus", "Rakiura", "Scirtidae", "Spaniocerca", "Spaniocercoides", "Zelandobius", "Zelolessica", "Zephlebia")
# write taxa forbid as table for MS
data.frame(Taxa = taxa.forbid) %>%
  write.csv("figures/Niche_forbid_taxa.csv",
            row.names = FALSE)

# set forbidden taxa columns = 0
P.niche <- map(P.niche,
         rm_niche,
         taxa = taxa.forbid)

# P.niche * relative abundance matrices
Pij <- map2(P.niche, N, ~.x*.y)

# rescale probabilities between 0 and 1
Pij <- map(Pij, scalexy, min = 0.01, max = 0.99)

# map(Pij, plot_heat)
# #map(Pij, hist)
# map(Pij, mean)

# save final probabilities which account for
# niche, neutral, and forbidden links
saveRDS(Pij, "data/AMD_final_probability_matrices.RDS")
