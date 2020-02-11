# estimate niche link probability
# June 2018
# jfpomeranz@gmail.com

# infer probability of links occuring in AMD dataset 
# based on body size, and parameters calculated using the integrated model
# script 2_trait_match_model.R

library(plyr)
library(dplyr)
library(traitmatch)
library(gplots)
library(purrr)
# function from Bartomeus
source("functions/predict.niche.prob.R")
#functions written for this ms
source("functions/MS_functions.R")

# model paramaters calculated using individual level data from Woodward et al. 2010. 
# for access to data for Broadstone Stream and Tadnoll Brook, please contact Guy Woodward and Iwan Jones, respectively
pars_pre <- readRDS("data/Integrated_model_params.RDS")

# AMD data ####
dataset <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
# arrange data by site and increasing body size
dataset <- arrange(dataset, site, avg.dw)
# split into list by site
dataset <- split(dataset, list(dataset$site))
# create df with all possible pairwise combinations of
# body sizes at each site in AMD data
M <- llply(dataset, function (x){
  expand.grid(log10(x$avg.dw), log10(x$avg.dw))
})

# estimate link probabilities
link.probs <- llply(M, function (x){
  predict.niche.prob(pars_pre,
                     x[[1]],
                     x[[2]],
                     replicates = 1)
})

# plot link probabilities for illustration
mybreaks <- seq(0,1,0.05)
for (i in 1:length(link.probs)){
  heatmap.2(matrix(link.probs[[i]][[1]],
                   length(dataset[[i]]$avg.dw),
                   length(dataset[[i]]$avg.dw)),
            Rowv = NA,
            Colv = NA,
            scale = "none",
            trace = "none",
            dendrogram = "none",
            breaks = mybreaks,
            key = F,
            labRow = NA,
            labCol = NA,
            main= c("Probability of interaction for site:",
                    names(link.probs)[i]),
            xlab = "Consumer",
            ylab = "Resource")
}

# probability matrix ####
# link.probs[[i]][[1]] == vector of link probabilities
# pull out just that element
prob.vec <- lapply(link.probs, "[[", 1)

# convert vector to a matrix. Vector is organized by dw
# as long as use the right ncol/nrow argument, matrix will be size sorted (e.g. typical predation matrix)

# convert all vectors to matrices
prob.matr <- map(prob.vec, get_prob_matr)
# add dimnames to matrices
for(i in 1:length(prob.matr)){
  dimnames(prob.matr[[i]]) <- list(dataset[[i]]$taxa,
                                   dataset[[i]]$taxa)
}
saveRDS(prob.matr, "data/AMD_niche_link_probability_matrices.RDS")
