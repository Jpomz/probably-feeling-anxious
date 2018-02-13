# 4_J_x_N

# estimate link probability by multiplying niche probability (J) and neutral probability (N)
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

# data ####
# read in link probabilities (J)
J <- readRDS("data/AMD link probability matrices.rds")
# read in relative abundance matrices
N <- readRDS("data/AMD relative abundance matrices.rds")

# scale N first
# scale N [0.5,1]
N.scale <- map(N, scalexy, min = 0.6, max = 1)
Aij <- map2(J, N.scale, ~.x*.y)
map(Aij, plot_heat)
map(Aij, hist)

map(Aij, mean)

b_trial <- function (prob_matr, trials){
  prob <- as.vector(prob_matr)
  out <- NULL
  for (i in 1:trials){
    out[[i]] <- matrix(rbinom(length(prob),
                       1, prob = prob),
                       ncol = ncol(prob_matr),
                       nrow = nrow(prob_matr))
  }
  out
}

itali <- b_trial(Aij$Italia, 100)
plot_heat(Aij$Italia)
points(itali[[1]])
