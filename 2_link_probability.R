# 1_link_probability
# 9 feb 2018
#J Pomz
# jfpomeranz@gmail.com

library(plyr)
library(tidyverse)
library(traitmatch)
# function from Bartomeus
source("predict.niche.prob.R")

bs <- read_csv("data/raw data/BroadstoneData_AERv45ch3.csv",
               skip = 16)[,-1]
# is bs in g or mg??? ####
# /1000
bs.g <- bs %>% mutate(predMass = predMass / 1000, preyMass = preyMass /1000)

war.all <- read_csv("data/raw data/warburton.csv")
# need to filter out observations when no predation occured
# e.g. "number_eaten_hr" == 0
# estimate dry mass
# (base^(ln_a + (b * log(length)))) / g_conversion
war.g <- war.all %>% 
  filter(number_eaten_hr >0) %>%
  transmute(predMass = (pred_base^
  (pred_ln_a + (pred_b * (log(pred_length_mm,
            base = pred_base))))) / pred_g_conversion,
            preyMass = prey_base^
    (prey_ln_a + (prey_b * (log(prey_length_mm,
            base = prey_base)))) / prey_g_conversion)

taieri <- readRDS("data/Taieri pred prey estimated dw g.RDS")
taieri <- rename(taieri, predMass=pred, preyMass=prey)
taieri$study <- "thompson"

all.dat <- bind_rows(war.g, bs.g, taieri)

ggplot(all.dat, aes(log10(predMass), log10(preyMass), color = study))+
  geom_point()+
  stat_quantile(quantiles = c(0.01,0.97))

# Trait match####
MPred <- log10(all.dat$predMass)
#c(log10(dat$predMass), log10(t.con*1000))
MPrey <- log10(all.dat$preyMass) 
#c(log10(dat$preyMass), log10(t.res*1000)) 

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
pars_pre
plot_pred(pars = pars_pre, 
          Tlevel1 = MPrey, 
          Tlevel2 = MPred, 
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)", 
          pch = "0")

pars_niche <- fit_it(niche_model, 
              Tlevel1 = MPrey,  
              Tlevel2 = MPred,
              mean_Tlevel1 = mean(MPrey),
              sd_Tlevel1 = sd(MPrey),
              pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
              par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
              par_hi = c(a0 = 10, a1 = 10, a1b0 = 10, b1 = 10),
              max.time = mt)
pars_pre
pars_niche
plot_pred(pars = pars_niche,
          Tlevel1 = MPrey, 
          Tlevel2 = MPred,
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)",
          pch = "0")

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

