# 1_link_probability
# 9 feb 2018
#J Pomz
# jfpomeranz@gmail.com

# 6 jun 2018
# using individual data from Broadstone stream and Tadnoll Brook
# received via email from Guy Woodward and Iwan Jones respectively

library(plyr)
library(tidyverse)
library(traitmatch)
library(gplots)
# function from Bartomeus
source("predict.niche.prob.R")

# broadstone ####
# individual pred-prey biomass data in mg
bs <- read_csv("data/raw data/BroadstoneData_AERv45ch3.csv", skip = 16)[,-1]
# convert to grams to match other datasets
bs.g <- bs %>% mutate(predMass = predMass / 1000,
                      preyMass = preyMass /1000)
bs.g$study <- "Woodward"

# tadnoll brook data in mg dw
tb <- read_csv("data/raw data/tadnoll pred prey biomass.csv")
#tb <- read_csv("data/raw data/TB_all_events.csv")[,-3]
tb <- tb[complete.cases(tb),]
# convert to grams
tb.g <- tb %>% mutate(predMass = predMass / 1000,
                      preyMass = preyMass /1000)
tb.g$study <- "tadnoll"

# combine data sets
all.dat <- bind_rows(bs.g, tb.g)
ggplot(all.dat, aes(log10(predMass), log10(preyMass), color = study))+
  geom_point()+
  stat_quantile(quantiles = c(0.05,0.95)) +
  theme_bw()

# Trait match####
MPred <- log10(all.dat$predMass)
MPrey <- log10(all.dat$preyMass) 

mt <- 60 #Define max.time to 60 sec to run things fast. Set to minimum 900 for a decent estimation of parameters.
pars_pre <- fit_it(integrated_model, 
              Tlevel1 = MPrey,  
              Tlevel2 = MPred,
              mean_Tlevel1 = mean(MPrey),
              sd_Tlevel1 = sd(MPrey),
              pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
              par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
              par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
              max.time = mt)
#pars_pre
plot_pred(pars = pars_pre,
          Tlevel1 = MPrey,
          Tlevel2 = MPred,
          xlab = "log (Predator body size)",
          ylab = "log (prey body size)",
          pch = "0")
# 
pars_niche <- fit_it(niche_model,
              Tlevel1 = MPrey,
              Tlevel2 = MPred,
              mean_Tlevel1 = mean(MPrey),
              sd_Tlevel1 = sd(MPrey),
              pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
              par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
              par_hi = c(a0 = 10, a1 = 10, a1b0 = 10, b1 = 10),
              max.time = mt)
# pars_pre
# pars_niche
plot_pred(pars = pars_niche,
          Tlevel1 = MPrey,
          Tlevel2 = MPred,
          xlab = "log (Predator body size)",
          ylab = "log (prey body size)",
          pch = "0")

# compare integrated, niche, neutral models
lh_model <- -integrated_model(pars_pre, MPrey,
                              MPred, mean(MPrey), sd(MPrey))
lh_niche <- -niche_model(pars_niche, MPrey, MPred,
                         mean(MPrey), sd(MPrey))
lh_neutral <- -neutral_model(pars = NULL, MPrey, MPred,
                             mean(MPrey), sd(MPrey))
barplot(c(lh_model, lh_niche, lh_neutral),
        names.arg = c("integrated",
                      "niche",
                      "neutral"))
# integrated is "better"

# AMD data ####
dataset <- readRDS("data/AMD fish invert dw abundance.RDS")
dataset <- arrange(dataset, site, avg.dw)
dataset <- split(dataset, list(dataset$site))
M <- llply(dataset, function (x){
  expand.grid(log10(x$avg.dw), log10(x$avg.dw))
})

link.probs <- llply(M, function (x){
  predict.niche.prob(pars_pre,
                     x[[1]],
                     x[[2]],
                     replicates = 100)
})

# plot link probabilities for illustration
mybreaks <- seq(0,1,0.1)
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
            main= "Probability of interaction",
            xlab = "Consumer",
            ylab = "Resource")
}

# link.probs[[i]][[1]] == vector of link probabilities
# pull out just that element
prob.vec <- lapply(link.probs, "[[", 1)
# convert vector to a matrix. Vector is organized by dw
# as long as use the right ncol/nrow argument, matrix will be size sorted (e.g. typical predation matrix)
# function to convert to matrix
prob_matr <- function(vec){
  m = matrix(vec, sqrt(length(vec)), sqrt(length(vec)))
}
# convert all vectors to matrices
prob.matr <- map(prob.vec, prob_matr)
# add dimnames to matrices
for(i in 1:length(prob.matr)){
  dimnames(prob.matr[[i]]) <- list(dataset[[i]]$taxa,
                                   dataset[[i]]$taxa)
}
saveRDS(prob.matr, "data/AMD link probability matrices.RDS")
