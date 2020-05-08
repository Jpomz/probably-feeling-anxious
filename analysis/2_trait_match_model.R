# 2_trait_match_model
# 9 feb 2018
# J Pomz on GitHub
# jfpomeranz@gmail.com

# 6 jun 2018
# using individual data from Broadstone stream and Tadnoll Brook
# received via email from Guy Woodward and Iwan Jones respectively
# please contact them directly if you would like to use this data

library(plyr)
library(tidyverse)
# to install "traitmatch" package, see:
# https://github.com/ibartomeus/traitmatch
# note that the syntax of the command:
# install_github("traitmatch", "ibartomeus")
# is deprecated. 
# replace with :
# install_github("ibartomeus/traitmatch")

library(traitmatch)

# function to predict new interaction probabilities 
# received from Ignasi Bartomeus
# if using, please cite traitmatch package above and Bartomeus et al. 2016. A common framework for identifying linkage rules across different types of interactions. Functional Ecology
source("functions/predict.niche.prob.R")


# broadstone ####
# individual pred-prey biomass data in mg
bs <- read_csv("data/raw data/BroadstoneData_AERv45ch3.csv", skip = 16)[,-1]
bs <- bs %>%
  select(predMass, preyMass)
# tadnoll brook data
# individual pred-prey biomass data in mg
tb <- read_csv("data/raw data/tadnoll pred prey biomass.csv")
# combine data sets
all.dat <- bind_rows(bs, tb)

# convert to grams to match other datasets
all.dat <- all.dat %>%
  mutate(predMass = predMass / 1000,
         preyMass = preyMass /1000)


ggplot(all.dat, aes(log10(predMass), log10(preyMass)))+
  geom_point(alpha = 0.2)+
  stat_quantile(quantiles = c(0.025,0.975)) +
  theme_bw()

# Trait match####
# see Bartomeus et al. 2015 for discussion on "integrated model"
MPred <- log10(all.dat$predMass)
MPrey <- log10(all.dat$preyMass) 

mt <- 900 #Define max.time to 60 sec to run things fast. Set to minimum 900 for a decent estimation of parameters.
# estimate parameters using integrated model
pars_integrated <- fit_it(integrated_model, 
              Tlevel1 = MPrey,  
              Tlevel2 = MPred,
              mean_Tlevel1 = mean(MPrey),
              sd_Tlevel1 = sd(MPrey),
              pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
              par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
              par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
              max.time = mt)
# plot integrated model
plot_pred(pars = pars_integrated,
          Tlevel1 = MPrey,
          Tlevel2 = MPred,
          xlab = "log (Predator body size)",
          ylab = "log (prey body size)",
          pch = "0")

saveRDS(pars_integrated, "data/Integrated_model_params.RDS")


# compare to other model fits
pars_pre_niche <- fit_it(niche_model,
                         Tlevel1 = MPrey, 
                         Tlevel2 = MPred,
                         mean_Tlevel1 = mean(MPrey),
                         sd_Tlevel1 = sd(MPrey),
                         pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                         par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
                         par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
                         max.time = mt)

lh_model <- -integrated_model(pars_integrated, MPrey, MPred, mean(MPrey),sd(MPrey))
lh_niche <- -niche_model(pars_pre_niche, MPrey, MPred, mean(MPrey), sd(MPrey))
lh_neutral <- -neutral_model(pars = NULL, MPrey, MPred, mean(MPrey), sd(MPrey))

# Visualization
barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))



# percent of realized interactions greater than some probability
# Calculate interaction probabilies (pLM) using function predict.niceh.prob()
int_prob <- predict.niche.prob(pars_pre, MPrey, MPred, replicates = 1)[[2]]

# column "pLM" is the predicted probability of an interaction from the model
int_prob %>%
  filter(pLM <0.5) %>% #probabilities > or < some probability e.g. 0.5
  count %>% # no. obs for given probabilty
  mutate(n / nrow(niche_probs)) # proportion of realized links for given probability

