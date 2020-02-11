# bayesian stability results
library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms)

theme_set(theme_bw())

#read in simulation data and split by .id
dat.df <- readRDS("data/stability_results.RDS")
# only keeping some variables for Ram space
dat.df <- dat.df %>% 
  select(Site, pca1, stab, .id)
dat <- split(dat.df, list(dat.df$.id))
# can also rm() dat.df to save more
rm(dat.df)

# read in model results
rand.brm <- readRDS("data/Bayesian_stability_random.RDS")
scale.brm <- readRDS("data/Bayesian_stability_scaled.RDS")
corr.brm <- readRDS("data/Bayesian_stability_corr.RDS")
s.c.brm <- readRDS( "data/Bayesian_stability_s_c.RDS")

# print(rand.brm)
# print(scale.brm)
# print(corr.brm)
# print(s.c.brm)

# plot model fits with raw data as points
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

plot_marg_eff(model = rand.brm,
              raw_data = dat$random,
              title = "Random interaction strengths") 
plot_marg_eff(model = scale.brm,
              raw_data = dat$scaled,
              title = "Scaled interaction strengths")
plot_marg_eff(model = corr.brm,
              raw_data = dat$corr,
              title = "Correlated interaction strengths")
plot_marg_eff(model = s.c.brm,
              raw_data = dat$s.c,
              title = "Scaled & Correlated interaction strengths")


# probability slope !=0
# using posterior_samples() to take total number of obs that != 0, dividing by total number of observations
# Is this what you did in your 2017 ES&T paper?
# not sure if this is necessary, but just trying to learn more by parroting your 2017 paper
coef_prob_not_0 <- function(model, coeficient){
  post_samples <- posterior_samples(model)[,coeficient]
  probability <- sum(post_samples!=0) / length(post_samples)
  return(probability)
}
coef_prob_not_0(rand.brm, "b_pca1")
  # = 1
  # erego, probability that slope of random model !=0 is 1???
coef_prob_not_0(rand.brm, "b_Intercept")
coef_prob_not_0(scale.brm, "b_pca1")
coef_prob_not_0(scale.brm, "b_Intercept")
coef_prob_not_0(corr.brm, "b_pca1")
coef_prob_not_0(corr.brm, "b_Intercept")
coef_prob_not_0(s.c.brm, "b_pca1")
coef_prob_not_0(s.c.brm, "b_Intercept")
### All of these are 1... so either good models or code is wrong?

plot_posterior <- function (model, raw_data, title = NULL, rel_min_height = 0){
  newdata <- raw_data %>% distinct(pca1, Site)
  fitted_new <- fitted(model,
                       newdata = newdata,
                       re_formula = NA,
                       summary = FALSE)
  colnames(fitted_new) <-  newdata %>% 
    pull(Site) #rename columns as the sites that correspond to pca1
  fitted_plot <- as.data.frame(fitted_new)
  
  fitted_plot %>% 
    gather(Site, stability) %>%
    merge(newdata) %>% 
    ggplot(aes(x = stability,
               y = reorder(Site, pca1),
               fill = reorder(Site, pca1)))+
    geom_density_ridges(rel_min_height = rel_min_height)+
    ggtitle(paste("Posterior distribution of stability as predicted from\n", title, "model"))+
    scale_fill_manual(values = colorRampPalette(rev(brewer.pal(
      n = 11,
      name = "Spectral")))(25)) +
    theme(legend.position = "NULL") +
    ylab("Site") +
    xlab(expression(paste("Stability, ", italic("s"))))
}

plot_posterior(rand.brm,
               raw_data = dat$random,
               title = "random interaction strength")
plot_posterior(scale.brm,
               raw_data = dat$scaled,
               title = "scaled interaction strength")
plot_posterior(corr.brm,
               raw_data = dat$corr,
               title = "correlated interaction strength")
plot_posterior(s.c.brm,
               raw_data = dat$s.c,
               title = "scaled & correlated interaction strength")


ggplot(dat.df, aes(x = stab,
                y = reorder(Site,pca1),
                fill = reorder(Site,pca1))) +
  geom_density_ridges(rel_min_height = 0.05) +
  # rel_min_height filters out bottom x.xx% 
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(
    n = 11,
    name = "Spectral")))(25)) +
  theme_bw() +
  theme(legend.position = "NULL") +
  labs(title = "Simulated stability metric",
       y = "Site",
       x = expression(italic("s"))) +
  facet_wrap(~.id)+
  coord_cartesian(xlim = c(0, 0.11))

