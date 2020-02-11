# interpret bayesian model outputs

# bayesian stability results
library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms)

theme_set(theme_bw())

# bayesian model outputs 
rand.brm <- readRDS("data/Bayesian_stability_random.RDS")
# scale.brm <- readRDS("data/Bayesian_stability_scaled.RDS")
# corr.brm <- readRDS("data/Bayesian_stability_corr.RDS")
# s.c.brm <- readRDS( "data/Bayesian_stability_s_c.RDS")
#============================
# marginal efffects plot ####
#============================

# simulation output data
# needed for points on plot
dat.df <- readRDS("data/stability_results.RDS")
# only keeping some variables for Ram space
dat.df <- dat.df %>% 
  select(Site, pca1, stab, .id)
dat <- split(dat.df, list(dat.df$.id))

#==================================
# plot marg eff for reference ####
#==================================
# function for plot
plot_marg_eff <- function(model, raw_data, title = NULL){
  marg_eff <- marginal_effects(model, method = "fitted")
  marg_eff_df <- as.data.frame(marg_eff$pca1)
  #return(marg_eff_df)
  ggplot()+
    geom_point(data = raw_data,
               aes(x = pca1,
                   y = stab,
                   color = reorder(Site, -pca1)
               ),
               alpha = 0.3,
               position = position_jitter(width=0.1))+
    geom_ribbon(data = marg_eff_df,
                aes(x = pca1,
                    ymin = lower__,
                    ymax = upper__),
                alpha = 0.15)+
    geom_line(data = marg_eff_df,
              aes(x = pca1,
                  y = estimate__),
              color = "blue",
              size=1) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),
                                 title="Site"))+
    scale_color_manual(values = colorRampPalette((brewer.pal(
      n = 11,
      name = "Spectral")))(25))+
    ylab(expression(paste("Stability, ", italic("s")))) +
    xlab("PCA1 (AMD Stress)") +
    ggtitle(title)
}

# plot stability values based on random interaction strengths
plot_marg_eff(model = rand.brm,
              raw_data = dat$random,
              title = "Random interaction strengths") 

#==================================
# random interaction strengths ####
#==================================

# check posterior predictive
pp_check(rand.brm, type="boxplot")
#looks good to JPZ

# extract posterior samples
post_random <- posterior_samples(rand.brm)

# mean 95% CrI for parameter values
# This should be exponentiated correct??
# e.g. NOT quantile(post_random$b_Intercept, probs=c(0.025,0.5,0.975))

# intercept
quantile(exp(post_random$b_Intercept), probs=c(0.025,0.5,0.975))
# the intercept value here (0.028) makes sense based on marginal effects plot (below)

# slope
quantile(exp(post_random$b_pca1), probs=c(0.025,0.5,0.975))
# I'm confused about this, shouldn't the slope be negative? Is the code above correct??

# confused about these too, see note above
# probability that slope is not 0
sum(exp(post_random$b_pca1) !=0) / nrow(post_random)
# = 0.9999999 

# probability that slope is < 0
sum(exp(post_random$b_pca1) < 0) / nrow(post_random)
# = 0

# probability that slope is > 0
sum(exp(post_random$b_pca1) >0) / nrow(post_random)
# = 0.9999
# again, I thought slope should be 0, so these probabilities seem backwards to me unless I'm missing something?
# see plot of marginal effects below

#============================
# relative changes across gradient ####
#============================

# testing three spots
# 1) reference to min_impacts (where AMD "starts")
# 2) reference to max_impacts
# 3) min_impacts to max_impacts

# pca1 values for three end poinst
# reference = -3.5599
# min_impacts = -0.8089
# max_impacts = 6.1127
reference = -3.5599
min_impacts = -0.8089
max_impacts = 6.1127
# calculate stability at 3 points
ref_stab <- exp(post_random$b_Intercept + post_random$b_pca1*reference)
min_impacts_stab <- exp(post_random$b_Intercept + post_random$b_pca1*min_impacts)
max_impacts_stab <- exp(post_random$b_Intercept + post_random$b_pca1*max_impacts)

# fold change with increasing impacts
# ref to min_impacts
quantile(ref_stab/min_impacts_stab, probs = c(0.025, 0.5, 0.975))
# stability values in reference sites are ~ 2x greater than sites with minimal impact (remember greater values indicate LOWER stability)

# ref to max_impacts
quantile(ref_stab/max_impacts_stab, probs = c(0.025, 0.5, 0.975))
# stability values in reference sites are ~ 14x greater than in sites with greatest impacts

# min to max impacts
quantile(min_impacts_stab/max_impacts_stab, probs = c(0.025, 0.5, 0.975))
# stability values in minimally impacted sites are 6x greater than in sites with greatest impacts

# wondering if these three calcs above should be reversed

