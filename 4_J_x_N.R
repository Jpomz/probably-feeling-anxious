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
            breaks = seq(0,1,0.1),
            key = F,
            labRow = NA,
            labCol = NA,
            main= "Probability of interaction",
            xlab = "Consumer",
            ylab = "Resource")
}

# function to rescale [0,1]
scale01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}
# function to rescale variable to [x,y]
scalexy <- function(x, min, max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

# data ####
# read in link probabilities (J)
J <- readRDS("data/AMD link probability matrices.rds")
# read in relative abundance matrices
N <- readRDS("data/AMD relative abundance matrices.rds")

# Aij ####
Aij <- map2(J, N, ~.x*.y)
map(Aij, plot_heat)
# scale Aij [0,1]
Aij.scale <- map(Aij, scale01)
map(Aij.scale, plot_heat)
# scale N first
N.scale <- map(N, scale01)
# aij with scaled N
Aij.n.scale <- map2(J, N.scale, ~.x*.y)

# scale N [0.5,1]
N.scalexy <- map(N, scalexy, min = 0.5, max = 1)
Aij.nxy <- map2(J, N.scalexy, ~.x*.y)
map(Aij.nxy, plot_heat)

map(J, plot_heat)
