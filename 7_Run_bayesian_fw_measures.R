# run Bayesian FW measures

library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms) 

# NOTE ================================#
# these analyses take a long time as well as a lot of memory
# You may need to restart your R session e.g. ctrl + shift+ F10 (on a PC) between each bayesian analysis. 

theme_set(theme_bw())
dat <- readRDS("data/stability_results.RDS")
dat <- dat %>% filter(.id == "random")

# sum(dat$L==0)
# sum(dat$C==0)


# sum(is.na(dat$mean.TL))
# sum(na.omit(dat$mean.TL)==0)
# sum(is.na(dat$max.TL))
# sum(na.omit(dat$max.TL)==0)
# sum(is.na(dat$sd.TL))
# sum(na.omit(dat$sd.TL)==0)

sum(dat$S==0)
dim(dat)
# #
S.brm <- brm(S ~ pca1 + (1|Site),
               data = dat,
               family = Gamma(link = "log"),
               prior=c(prior(normal(0, 1),class = "b"),
                       prior(normal(0, 1),class = "Intercept")),
               chains = 4,
               iter = 2000,
               cores = 4,
             control = list(max_treedepth = 15))
saveRDS(S.brm, "data/Bayesian_FW_S.RDS")
pairs(S.brm)


# L.brm <- brm(L ~ pca1 + (1|Site),
#                 data = dat,
#                 family = Gamma(link = "log"),
#                 prior=c(prior(normal(0, 1),class = "b"),
#                         prior(normal(0, 1),class = "Intercept")),
#                 chains = 4,
#                 iter = 2000,
#                 cores = 4)
# saveRDS(L.brm, "data/Bayesian_FW_L.RDS")

# C.brm <- brm(C ~ pca1 + (1|Site),
#              data = dat,
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(C.brm, "data/Bayesian_FW_C.RDS")

# sum(dat$Gensd==0)
# dim(dat)
# dat <- dat[-which(dat$Gensd==0),]
# sum(dat$Gensd==0)
# dim(dat)
# #

#Gen SD bayesian
# Gen.brm <- brm(Gensd ~ pca1 + (1|Site),
#              data = dat,
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(Gen.brm, "data/Bayesian_FW_Gen.RDS")
# pairs(Gen.brm)

# Calculate generality and vulnerability ####
# generality: G = L /(Nt + Ni)
# vulnerability: V = L / (Nb + Ni)
# L = Links, Nb, Ni, Nt = number of bottom, intermediate and top species, respectively. 
dat <- dat %>% 
  mutate(G = (L / (S* T+ S * I)) / S,
         V = (L / (S * B + S * I)) / S)

sum(dat$G == 0)

Gen.brm <- brm(G ~ pca1 + (1|Site),
             data = dat,
             family = Gamma(link = "log"),
             prior=c(prior(normal(0, 1),class = "b"),
                     prior(normal(0, 1),class = "Intercept")),
             chains = 4,
             iter = 2000,
             cores = 4)
saveRDS(Gen.brm, "data/Bayesian_FW_Gen.RDS")
pairs(Gen.brm)


sum(dat$V == 0)
Vul.brm <- brm(V ~ pca1 + (1|Site),
               data = dat,
               family = Gamma(link = "log"),
               prior=c(prior(normal(0, 1),class = "b"),
                       prior(normal(0, 1),class = "Intercept")),
               chains = 4,
               iter = 2000,
               cores = 4)
saveRDS(Vul.brm, "data/Bayesian_FW_Vul.RDS")
pairs(Vul.brm)

# sum(dat$Vulsd==0)
# 
# Vul.brm <- brm(Vulsd ~ pca1 + (1|Site),
#                data = dat,
#                family = Gamma(link = "log"),
#                prior=c(prior(normal(0, 1),class = "b"),
#                        prior(normal(0, 1),class = "Intercept")),
#                chains = 4,
#                iter = 2000,
#                cores = 4)
# saveRDS(Vul.brm, "data/Bayesian_FW_Vul.RDS")

# B.brm <- brm(B ~ pca1 + (1|Site),
#              data = dat,
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(B.brm, "data/Bayesian_FW_B.RDS")
# 
# I.brm <- brm(I ~ pca1 + (1|Site),
#              data = dat,
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(I.brm, "data/Bayesian_FW_I.RDS")
# 
# T.brm <- brm(T ~ pca1 + (1|Site),
#              data = dat,
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(T.brm, "data/Bayesian_FW_T.RDS")

# sum(dat$max.TL==0)
# 
# max.TL.brm <- brm(max.TL ~ pca1 + (1|Site),
#              data = na.omit(dat),
#              family = Gamma(link = "log"),
#              prior=c(prior(normal(0, 1),class = "b"),
#                      prior(normal(0, 1),class = "Intercept")),
#              chains = 4,
#              iter = 2000,
#              cores = 4)
# saveRDS(max.TL.brm, "data/Bayesian_FW_max.TL.RDS")
# 
# mean.TL.brm <- brm(mean.TL ~ pca1 + (1|Site),
#                   data = na.omit(dat),
#                   family = Gamma(link = "log"),
#                   prior=c(prior(normal(0, 1),class = "b"),
#                           prior(normal(0, 1),class = "Intercept")),
#                   chains = 4,
#                   iter = 2000,
#                   cores = 4)
# saveRDS(mean.TL.brm, "data/Bayesian_FW_mean.TL.RDS")
