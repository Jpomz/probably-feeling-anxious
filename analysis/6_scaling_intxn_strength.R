# Different interaction strength estimates
# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(forcats)
library(gplots)

# food web functions from Petchey
source("functions/FoodWebFunctions.r")
# functions written for this MS
source("functions/MS_functions.R")
# stability functions from Sauve 2016
# if using, please cite original publication:
# Sauve, A. M. C., Thébault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1
source("functions/stability_fxns_Sauve.R")


# read in data
# probabilty matrices
Pij <- readRDS("data/AMD_final_probability_matrices.RDS")
# gradient
gradient <- readRDS("data/pca_axis_26_sites.RDS")


# random interaction strengths ####
# set seed for reproducibility
set.seed(3049)
random.J <- ldply(
  map(
    map(
      Pij,
      get_measures, 
      s2 = 2,
      trials = 250),
    ldply)
)
# join gradient + random.J
random.J <- left_join(random.J, gradient, by = c(".id" = "site"))

# scaled interaction strengths ####
# for reproducibility
set.seed(3049)
scaled.J <- ldply(
  map(
    map(
      Pij,
      get_measures,
      s2 = 2,
      trials = 250,
      scale.Jij = TRUE),
    ldply)
  )

scaled.J <- left_join(scaled.J, gradient, by = c(".id" = "site"))

# correlated interaction strengths ####
# for reproducibility
set.seed(3049)
corr.J <- ldply(
  map(
    map(
      Pij,
      get_measures, 
      s2 = 2, 
      trials = 250, 
      correlate.Jij = TRUE),
    ldply)
)
corr.J <- left_join(corr.J, gradient, by = c(".id" = "site"))

# scaled and correlated interaction strengths ####
# for reproducibility
set.seed(3049)
s.c.J <- ldply(
  map(
    map(
      Pij,
      get_measures,
      s2 = 2,
      trials = 250, 
      correlate.Jij = TRUE,
      scale.Jij = TRUE),
    ldply),
)

s.c.J <- left_join(s.c.J, gradient, by = c(".id" = "site"))

# combine all stability results into a list. 
# removing column 1 (".id") from each dataframe. This information is redundant with "Site" column
stability.results <- list(random = random.J[,-1],
                          scaled = scaled.J[,-1],
                          corr = corr.J[,-1],
                          s.c = s.c.J[,-1])
stability.results <- ldply(stability.results)

saveRDS(stability.results, "data/stability_results.RDS")
