# Bayesian attempt
library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms) 

dat <- readRDS("data/stability_results.RDS") # will need to update file location with wherever you save .RDS file

# dat has 4 different stability estimates specified by $.id column

# filter out one type for example 
# (eventually I would like to run this analysis for all 4 types)
rand.j <- dat %>% filter(.id == "random")
# random = "random" interaction strengths
# stab = qualitative stability measure (lower = "more stable")
# pca1 = principle component axis 1 e.g. continuous measure of AMD impact; higher values = increasig impact
# Site = site name, factor
  # should site be an ordered factor?

min(rand.j$stab)
# stab contains 0's
# this is artifact of inferring stability

sum(rand.j$stab==0)
# there are 13 observations where stab = 0 (out of 6250 total observations)

rand.j$stab[rand.j$stab==0] <- quantile(rand.j$stab, probs = 0.01)
# replace with lowest values observed in other sites 
# (quantile .01 to ~.20 == 6.1e-05)
# could also remove these observations if that is more statistically appropriate
sum(rand.j$stab==0)
# 0's replaced

# joyplot of modeled stability metric
ggplot(rand.j, aes(x = stab,
                   y = fct_reorder(Site, pca1),
                   fill = fct_reorder(Site, pca1))) +
  geom_density_ridges(rel_min_height = 0.05) +
  # rel_min_height filters out bottom x.xx% 
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(
    n = 11,
    name = "Spectral")))(25)) +
  theme_bw() +
  theme(legend.position = "NULL") +
  labs(title = "Random Jij",
       y = "Site",
       x = "s") +
  coord_cartesian(xlim = c(0, 0.11))
# sites increas in AMD stress bottom to top
# stability increases R to L
# e.g. impacted sites "more" stable


# Bayesian model with pca1 as fixed effect and site as random effect [Is this correct syntax?]
rand.brm <- brm(stab ~ pca1 + (1|Site),
                 data = rand.j,
                 family = Gamma(link = "log"),
                 prior=c(prior(normal(0, 1),class = "b"),
                         prior(normal(0, 1),class = "Intercept")),
                 chains=2,
                 iter=1000,
                 cores=4)

print(rand.brm)
# interpret output:
# what does the group-level effects mean? 

# pop level:
# pca1 = -0.27 --> with every 1 unit increase in pca1 gradient, the stability metric decreases by 0.27 [CrI -0.43, -0.14] ??

# family specific:
# no idea what this means

summary(rand.brm) # looks to be same as print?

plot(rand.brm)
# are the left panels the posterior distributions of the parameter?

marginal_effects(rand.brm, method="fitted")
# this looks to be a quadratic response
# is there any way to test that? or does this suggest we should actually look at a non-linear brm() e.g. stab ~ pca1 + I(pca1^2) + (1|Site)

marg_eff <- marginal_effects(rand.brm, method = "fitted")

marg_eff$pca1 # data frame with predicted stability and CrI across pca1 gradient
# no Site level info?

# from here I get a bit lost and not sure where to go. 
# having joyplots of posterior dists would be cool (similar to figure at beginning of script), but not necessary if there is a better way. 


