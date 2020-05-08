# Bayesian Stability Stats
# Justin Pomeranz and Jeff Wesner
# load libraries
library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms) 

theme_set(theme_bw())
# load stability results
# dataframe with results for interaction strength that are random, scaled, etc. 
dat <- readRDS("data/stability_results.RDS")

# split df into a list based on interaction strength estimate
dat <- split(dat, list(dat$.id))

# check for 0's
sum(dat$random$stab==0)
sum(dat$corr$stab==0)
sum(dat$scaled$stab==0)
sum(dat$s.c$stab==0)

# remove 0's
dat$random <- dat$random[-which(dat$random$stab==0),]
dat$corr <- dat$corr[-which(dat$corr$stab==0),]
dat$scaled <- dat$scaled[-which(dat$scaled$stab==0),]
dat$s.c <- dat$s.c[-which(dat$s.c$stab==0),]

# check that 0's are gone
sum(dat$random$stab==0)
sum(dat$corr$stab==0)
sum(dat$scaled$stab==0)
sum(dat$s.c$stab==0)

# Random interaction strength model
rand.brm <- brm(stab ~ pca1 + (1|Site),
                data = dat$random,
                family = Gamma(link = "log"),
                prior=c(prior(normal(0, 1),class = "b"),
                        prior(normal(0, 1),class = "Intercept")),
                chains = 4,
                iter = 2000,
                cores = 4)
saveRDS(rand.brm, "data/Bayesian_stability_random.RDS")

# Scaled interaction strength model
scale.brm <- brm(stab ~ pca1 + (1|Site),
                data = dat$scaled,
                family = Gamma(link = "log"),
                prior=c(prior(normal(0, 1),class = "b"),
                        prior(normal(0, 1),class = "Intercept")),
                chains = 4,
                iter = 2000,
                cores = 4)
saveRDS(scale.brm, "data/Bayesian_stability_scaled.RDS")

# correlated interaction strength model
corr.brm <- brm(stab ~ pca1 + (1|Site),
                 data = dat$corr,
                 family = Gamma(link = "log"),
                 prior=c(prior(normal(0, 1),class = "b"),
                         prior(normal(0, 1),class = "Intercept")),
                 chains = 4,
                 iter = 2000,
                 cores = 4)
saveRDS(corr.brm, "data/Bayesian_stability_corr.RDS")

# Scaled and correlated interaction strength model
s.c.brm <- brm(stab ~ pca1 + (1|Site),
                 data = dat$s.c,
                 family = Gamma(link = "log"),
                 prior=c(prior(normal(0, 1),class = "b"),
                         prior(normal(0, 1),class = "Intercept")),
                 chains = 4,
                 iter = 2000,
                 cores = 4)
saveRDS(s.c.brm, "data/Bayesian_stability_s_c.RDS")

