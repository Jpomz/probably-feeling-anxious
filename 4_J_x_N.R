# 4_J_x_N

# estimate link probability by multiplying niche probability (J) and neutral probability (N)
# also prune niche forbidden links
# from Poisot et al 2014: 
# Aij [proportional to] J(i,j) * N(i,j)
# J = probability calculated by predict.niche.prob() from Bartomeus, in script #2_link_probability.R
# N ~ ni * nj where n is a vector of relative abundances. 
# N is proportional to ni*nj, and the actual values of N are very low. Need to rescale this variable to properly interpret. 

library(plyr)
library(tidyverse)
library(gplots)

# function to plot heat map 
# just using this for data exploration
plot_heat <- function(x, ...){
  heatmap.2(x,
            Rowv = NA,
            Colv = NA,
            scale = "none",
            trace = "none",
            dendrogram = "none",
            breaks = seq(0,1,0.01),
            key = F,
            labRow = NA,
            labCol = NA,
            main= "Probability of interaction",
            xlab = "Consumer",
            ylab = "Resource")
}


# function to rescale variable to [min, max]
scalexy <- function(x, min, max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

# function to set probabilities of forbidden taxa to 0
rm_niche <- function(matrix, taxa){
  for(name in (
    colnames(matrix)[colnames(matrix) %in% taxa])){
    matrix[,name] <- 0
  }
  matrix
}
# function to "fix" probabilities of fish predators
fish_prob <- function(matrix, taxa, mean, sd){
  for(name in (
    colnames(matrix)[colnames(matrix) %in% taxa])){
    matrix[,name] <- rnorm(nrow(matrix), mean, sd)
  }
  matrix
}

# data ####
# read in link probabilities (J)
J <- readRDS("data/AMD link probability matrices.rds")
# read in relative abundance matrices
N <- readRDS("data/AMD relative abundance matrices.rds")

# prune niche forbidden links
taxa.forbid <- c("Acari", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Copepoda", "Deleatidium", "Elmidae", "Elmidae Adult", "Eriopterini", "Helicopsyche", "Hexatomini", "Hydraenidae", "Hydrophilidae", "Nesameletus","Oligochaetae", "Olinga","Ostracoda", "Oxyethira", "Paraleptamphopus", "Platyhelminthes", "Potamopyrgus", "Rakiura", "Scirtidae", "Spaniocerca", "Spaniocercoides", "Zelandobius", "Zelolessica", "Zephlebia")
taxa.fish <- c("Salmo.trutta", "Anguilla.australis",
            "Galaxias.fasciatus", "Galaxias.maculatus",
            "Gobiomorphus.huttoni")
J <- map(J,
         rm_niche,
         taxa = taxa.forbid)
# fix fish probabilities rnorm(mean = 0.7, 0.1)
J <- map(J,
         fish_prob,
         taxa = taxa.fish,
         mean = 0.7, sd = 0.1)

# scale N first
# scale N [0.5,1]
N.scale <- map(N, scalexy, min = 0.5, max = 1)
# link probability = J x N
Aij <- map2(J, N.scale, ~.x*.y)

# rescale probabilities between 0 and 1
Aij <- map(Aij, scalexy, min = 0, max = 1)

#map(Aij, plot_heat)
#map(Aij, hist)
#map(Aij, mean)

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

# b.trials <- map(Aij, b_trial, trials = 100)
# 
# saveRDS(b.trials,
#         paste("data/100 bernouli trials ",
#               Sys.Date(),
#               ".rds",
#               sep = ""))


b.trials <- map(Aij, b_trial, trials = 500)
saveRDS(b.trials, "data/500 bernouli trials.rds")
