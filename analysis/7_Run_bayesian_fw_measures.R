# Justin Pomeranz and Jeff Wesner
# run Bayesian FW measures across AMD gradient

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

# Number of species 
S.dat <- dat %>%
  distinct(pca1, S, Site)

S.brm <- brm(S ~ pca1 + (1|Site),
               data = S.dat,
               family = Gamma(link = "log"),
               prior=c(prior(normal(0, 1),class = "b"),
                       prior(normal(0, 1),class = "Intercept")),
               chains = 4,
               iter = 2000,
               cores = 4,
             control = list(adapt_delta = 0.99))
saveRDS(S.brm, "data/Bayesian_FW_S.RDS")
pairs(S.brm)

# Number of links
L.brm <- brm(L ~ pca1 + (1|Site),
                data = dat,
                family = Gamma(link = "log"),
                prior=c(prior(normal(0, 1),class = "b"),
                        prior(normal(0, 1),class = "Intercept")),
                chains = 4,
                iter = 2000,
                cores = 4)
saveRDS(L.brm, "data/Bayesian_FW_L.RDS")

# Connectance values
C.brm <- brm(C ~ pca1 + (1|Site),
             data = dat,
             family = Gamma(link = "log"),
             prior=c(prior(normal(0, 1),class = "b"),
                     prior(normal(0, 1),class = "Intercept")),
             chains = 4,
             iter = 2000,
             cores = 4)
saveRDS(C.brm, "data/Bayesian_FW_C.RDS")

# Calculate generality and vulnerability ####
# generality: G = L /(Nt + Ni)
# vulnerability: V = L / (Nb + Ni)
# L = Links, Nb, Ni, Nt = number of bottom, intermediate and top species, respectively. 
dat <- dat %>% 
  mutate(G = (L / (S* T+ S * I)) / S,
         V = (L / (S * B + S * I)) / S)

# Bayesian model for G
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


# Model for V
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