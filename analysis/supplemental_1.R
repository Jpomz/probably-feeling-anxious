# Pomeranz, Wesner and Harding Ecology
# supplemental information
# Method for inferring food web structure

library(traitmatch)
library(gplots)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
# function from Bartomeus
source("functions/predict.niche.prob.R")
source("functions/MS_functions.R")
source("functions/FoodWebFunctions.r")


# set seed for reproducibility
set.seed(9677)

# load model parameters
# calculated in script 2
pars_pre <- readRDS("data/Integrated_model_params.RDS")

# read in Dempsters Adjacency matrix and taxa data
A <- as.matrix(read.csv("data/supplemental_data/dempsters_A.csv",
                        row.names = 1))
taxa <- as.vector(read.csv("data/supplemental_data/dempsters_taxa.csv",
                 stringsAsFactors = FALSE))
niche.forbid <- as.character(
  read.csv("data/supplemental_data/niche_forbidden_taxa.csv",
           stringsAsFactors = FALSE,
           header = FALSE)[[1]])

dw <- taxa$dw
ab <- taxa$no.m2 / sum(taxa$no.m2)

# expand dw vector to get all possible pairwise combinations
M <- expand.grid(log10(dw), log10(dw))
# calculate interaction probabilities for all species pairs
P <-  predict.niche.prob(pars_pre,
                         M[[1]], 
                         M[[2]],
                         replicates = 1)[[1]]

# convert probability vector to matrix
P <- get_prob_matr(P)
# add taxa names to P
colnames(P) <- colnames(A)
# example plot of probability matrix
plot_heat(P)

# remove niche forbidden taxa
P.niche <- rm_niche(P, taxa = niche.forbid)
plot_heat(P.niche)

# calculate relative abundance matrix
N <- get_rel_ab(ab, taxa$taxa)
# multiply fish abundance by correction factor (Pomeranz et al. 2019 Methods in Ecology and Evolution)
N <- correct_fish_rel_ab(N,
                         fish = c("Gobiomorphus", "Galaxias",
                                  "Salmo", "Anguilla"),
                         cf = 10^3)

# foodweb measures of empirical web
A_meas <- round(Get.web.stats(A, which.stats = 1), 5)
A_meas[10] <-(A_meas[2] / (A_meas[1] * A_meas[6] + A_meas[1] * A_meas[5])) / A_meas[1]
A_meas[11] <-(A_meas[2] / (A_meas[1] * A_meas[4] + A_meas[1] * A_meas[5])) / A_meas[1]

names(A_meas) <- c("S", "L", "C", "B", "I", "T", "U", "CB", "PCB", "G", "V")

# calculate density of inferred food web measures
get_inf_measures <- function(P, N, min.scale = 0.5, trials = 500){
  result <- list()
  N_scale <- scalexy(N, min = min.scale, max = 1)
  P_scale <- P * N_scale
  for (i in 1:trials){
    b <- b_trial(P_scale)
    result[[i]] <- round(Get.web.stats(b, which.stats = 1), 5)
  }
  result <- ldply(result)
  result
}

inf_meas_25 <- get_inf_measures(P = P, N = N, min.scale = 0.25)

ggplot(inf_meas_25, aes(L)) +
  geom_density() +
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Re-scale abundance 0.25 to 1")

inf_meas_40 <- get_inf_measures(P = P, N = N, min.scale = 0.4)

ggplot(inf_meas_40, aes(L)) +
  geom_density() +
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Re-scale abundance 0.4 to 1")

inf_meas_50 <- get_inf_measures(P = P, N = N, min.scale = 0.5)

ggplot(inf_meas_50, aes(L)) +
  geom_density() +
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "No. Links. N Rescaled from 0.5 to 1")

ggplot(inf_meas_50, aes(C)) +
  geom_density() +
  annotate("segment",
           x = A_meas[3], xend = A_meas[3],
           y = 5, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Connectance. N Rescaled from 0.5 to 1")

inf_meas_50 <- inf_meas_50 %>%
  # calculate normalized G and V
  mutate(G = (L / (S* T+ S * I)) / S,
         V = (L / (S * B + S * I)) / S)

ggplot(inf_meas_50, aes(G)) +
  geom_density() +
  annotate("segment",
           x = A_meas[10], xend = A_meas[10],
           y = 25, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Normalized Generality. N Rescaled to 0.5 to 1")

ggplot(inf_meas_50, aes(V)) +
  geom_density() +
  annotate("segment",
           x = A_meas[11], xend = A_meas[11],
           y = 25, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Normalized Vulnerability. N Rescaled to 0.5 to 1")


inf_meas_55 <- get_inf_measures(P = P, N = N, min.scale = 0.55)

ggplot(inf_meas_55, aes(L)) +
  geom_density() +
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Re-scale abundance 0.55 to 1")

inf_meas_75 <- get_inf_measures(P = P, N = N, min.scale = 0.75)
ggplot(inf_meas_75, aes(L)) +
  geom_density()+
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Re-scale abundance 0.75 to 1")

inf_meas_90 <- get_inf_measures(P = P, N = N, min.scale = 0.9)
ggplot(inf_meas_90, aes(L)) +
  geom_density()+
  annotate("segment",
           x = A_meas[2], xend = A_meas[2],
           y = 0.01, yend = 0,
           size = 1.5, colour = "black", arrow=arrow()) +
  ggtitle(label = "Re-scale abundance 0.9 to 1")

# plot final niche probability matrix
plot_heat(scalexy(P.niche * scalexy(N, min = 0.5, max = 1), min = 0.0, max = 1))

# function to calculate true skill statistic and its components
get_tss <- function (observed, inferred, all.measures = TRUE){
  # observed and inferred adjacency matrices
  # make sure adjacency matrices are same dimensions and have same colnames
  #stopifnot(dim(observed) == dim(inferred), 
  #identical(colnames(observed), colnames(inferred)))
  # subtract inferred from observed
  minus <- observed - inferred
  # multiply observed and inferred
  multiply <- observed * inferred
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  S2 = ncol(prediction)**2
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
  a = length(prediction[prediction==2])
  b = length(prediction[prediction==-1]) 
  c = length(prediction[prediction==1]) 
  d = length(prediction[prediction==0])
  # calculate TSS
  tss = (a*d - b*c)/((a+c)*(b+d))
  if (all.measures == TRUE){
    vars = data.frame(tp = a/S2,
                      fp = b / S2,
                      fn = c / S2,
                      tn = d / S2)
    tss = cbind(tss, vars)
    return(tss)
  } 
  tss
}



# function to measure tss and components of inferred food webs after 
# cumulative trials
cumulative_b_trials <- function(prob_matr, A_matr, N_matr, scale.min = 0.5){
  require(plyr)
  cum_b <- matrix(data = 0,
                  nrow = nrow(prob_matr),
                  ncol=ncol(prob_matr))
  false_neg <- 1
  trial <- 0
  dataout <- list()
  while (false_neg>0) {
    b <- b_trial(prob_matr*scalexy(N_matr, min = scale.min, max = 1))
    cum_b <- cum_b+b
    cum_b[cum_b>1] <- 1
    tss <- get_tss(A_matr, cum_b)
    trial <- trial+1
    vars <- cbind(tss, trial)
    dataout[[trial]] <- vars
    false_neg <- vars$fn
    if(trial == 1000){
      break
    }
  }
  ldply(dataout)
}

# investigate number of Bernoulli trials for inferring food web structure
tss_b_trials <- cumulative_b_trials(prob_matr = P,
                               A_matr = A,
                               N_matr = N,
                               scale.min = 0.5)


tss_b_trials <- gather(tss_b_trials, "variable","response", 1:5)

# plot cumulative TSS
tss_b_trials %>%
  filter(variable == "tss") %>%
  ggplot(aes(x = trial, y = response)) +
  geom_point()

ggplot(tss_b_trials, aes(x = trial, y = response, color = variable)) +
  geom_point(alpha = 0.2) +
  theme_bw()

tss_b_trials %>%
  filter(variable == "tp" | variable == "fn") %>%
  ggplot(aes(x = trial, y = response, color = variable)) +
  geom_point()

tss_b_trials %>%
  filter(variable != "tss" ) %>%
  ggplot(aes(x = trial, y = response, color = variable)) +
  geom_point()

tss_b_trials %>%
  filter(variable == "fp" | variable == "tn") %>%
  ggplot(aes(x = trial, y = response, color = variable)) +
  geom_point()
