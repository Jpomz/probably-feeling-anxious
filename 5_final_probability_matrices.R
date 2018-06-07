# 4_P.niche_x_N

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

# # function to plot heat map 
# # just using this for data exploration
# plot_heat <- function(x, ...){
#   heatmap.2(x,
#             Rowv = NA,
#             Colv = NA,
#             scale = "none",
#             trace = "none",
#             dendrogram = "none",
#             breaks = seq(0,1,0.01),
#             key = F,
#             labRow = NA,
#             labCol = NA,
#             main= "Probability of interaction",
#             xlab = "Consumer",
#             ylab = "Resource")
# }

# function to set probabilities of forbidden taxa to 0
rm_niche <- function(matrix, taxa){
  for(name in (colnames(matrix)[colnames(matrix) %in% taxa])){
    matrix[,name] <- 0
  }
  matrix
}

# function to rescale variable to [min, max]
scalexy <- function(x, min, max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

# data ####
# read in link probabilities (P.niche)
P.niche <- readRDS("data/AMD_niche_link_probability_matrices.RDS")
# read in relative abundance matrices
N <- readRDS("data/AMD_relative_abundance_matrices.RDS")

# niche forbidden ####
# "forbidden" taxa based on morphology
# e.g. "brushing/scraping" mouthparts, 
# specialized filtering, Austrosimulium's cephalic fan
taxa.forbid <- c("Acari", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Copepoda", "Deleatidium", "Elmidae", "Elmidae Adult", "Eriopterini", "Helicopsyche", "Hexatomini", "Hydraenidae", "Hydrophilidae", "Nesameletus","Oligochaetae", "Olinga","Ostracoda", "Oxyethira", "Paraleptamphopus", "Platyhelminthes", "Potamopyrgus", "Rakiura", "Scirtidae", "Spaniocerca", "Spaniocercoides", "Zelandobius", "Zelolessica", "Zephlebia")

# set forbidden taxa columns = 0
P.niche <- map(P.niche,
         rm_niche,
         taxa = taxa.forbid)

# P.niche * relative abundance matrices
Pij <- map2(P.niche, N, ~.x*.y)

# rescale probabilities between 0 and 1
Pij <- map(Pij, scalexy, min = 0.01, max = 0.99)

#map(Pij, plot_heat)
#map(Pij, hist)
#map(Pij, mean)

# save final probabilities which account for
# niche, neutral, and forbidden links
saveRDS(Pij, "data/AMD_final_probability_matrices.RDS")






b_trial <- function (prob_matr, trials){
  out <- NULL
  for (i in 1:trials){
    out[[i]] <- matrix(rbinom(length(prob_matr),
                       1, prob = prob_matr),
                       ncol = ncol(prob_matr),
                       nrow = nrow(prob_matr))
  }
  out
}

# b.trials <- map(Pij, b_trial, trials = 100)
# 
# saveRDS(b.trials,
#         paste("data/100 bernouli trials ",
#               Sys.Date(),
#               ".rds",
#               sep = ""))


b.trials <- map(Pij, b_trial, trials = 500)
saveRDS(b.trials, "data/500 bernouli trials.rds")
