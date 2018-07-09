# 2_trait_match_model
# 9 feb 2018
#J Pomz
# jfpomeranz@gmail.com

# 6 jun 2018
# using individual data from Broadstone stream and Tadnoll Brook
# received via email from Guy Woodward and Iwan Jones respectively

library(plyr)
library(tidyverse)
library(traitmatch)


# broadstone ####
# individual pred-prey biomass data in mg
bs <- read_csv("data/raw data/BroadstoneData_AERv45ch3.csv", skip = 16)[,-1]
bs <- bs %>%
  select(predMass, preyMass)
# tadnoll brook data in mg dw
tb <- read_csv("data/raw data/tadnoll pred prey biomass.csv")
# combine data sets
all.dat <- bind_rows(bs, tb)

# convert to grams to match other datasets
all.dat <- all.dat %>%
  mutate(predMass = predMass / 1000,
         preyMass = preyMass /1000)


ggplot(all.dat, aes(log10(predMass), log10(preyMass)))+
  geom_point()+
  stat_quantile(quantiles = c(0.05,0.95)) +
  theme_bw()

# Trait match####
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
# # plot integrated model
# plot_pred(pars = pars_integrated,
#           Tlevel1 = MPrey,
#           Tlevel2 = MPred,
#           xlab = "log (Predator body size)",
#           ylab = "log (prey body size)",
#           pch = "0")
# # estimate parameters using niche model
pars_niche <- fit_it(niche_model,
              Tlevel1 = MPrey,
              Tlevel2 = MPred,
              mean_Tlevel1 = mean(MPrey),
              sd_Tlevel1 = sd(MPrey),
              pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
              par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
              par_hi = c(a0 = 10, a1 = 10, a1b0 = 10, b1 = 10),
              max.time = mt)

# # Plot niche model
# plot_pred(pars = pars_niche,
#           Tlevel1 = MPrey,
#           Tlevel2 = MPred,
#           xlab = "log (Predator body size)",
#           ylab = "log (prey body size)",
#           pch = "0")

# compare integrated, niche, neutral models
lh_model <- -integrated_model(pars_integrated, MPrey,
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

saveRDS(pars_integrated, "data/Integrated_model_params.RDS")


