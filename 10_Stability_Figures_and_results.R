# interpret stability bayesian model outputs

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(brms)
source("functions/MS_functions.R")

theme_set(theme_bw())

# simulation output data
# needed for points on marginal effects plot
dat.df <- readRDS("data/stability_results.RDS")
# only keeping some variables for Ram space
dat.df <- dat.df %>% 
  select(Site, pca1, stab, .id)
dat <- split(dat.df, list(dat.df$.id))
# remove dat.df for ram space
rm(dat.df)


#============================
# stability ####
#============================

# Random interaction strengths ####
# read in Bayesian model for random interaction strengths
rand.brm <- readRDS("data/Bayesian_stability_random.RDS")

# marginal efffects plot ####
rand.marg.eff <- plot_marg_eff(model = rand.brm,
              raw_data = dat$random,
              title = "Random interaction strengths") 

# check posterior predictive ability of model
pp_check(rand.brm, type="boxplot")
# Boxplots largely overlap
# even large values, dots above boxplots are close
# this model can accurately predict observed values

# print out model summary
# (see MS_functions.R for more details)
rand.summary <- model_summary(rand.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/Random_model_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
    plot_prior_vs_posterior(rand.brm)
dev.off()

# remove randon model for RAM space
rm(rand.brm)

#============================
# scaled interaction strengths ####
scale.brm <- readRDS("data/Bayesian_stability_scaled.RDS")

# marginal efffects plot ####
scale.marg.eff <- plot_marg_eff(
  model = scale.brm,
  raw_data = dat$scaled,
  title = "Scaled interaction strengths") 

# check posterior predictive ability of model
pp_check(scale.brm, type="boxplot")
# this model also looks good
# print out model summary
# (see MS_functions.R for more details)
scale.summary <- model_summary(scale.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/Scaled_model_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
    plot_prior_vs_posterior(scale.brm)
dev.off()

# remove model for RAM space
rm(scale.brm)

#============================
# correlated interaction strengths ####
corr.brm <- readRDS("data/Bayesian_stability_corr.RDS")
# marginal efffects plot ####
corr.marg.eff <- plot_marg_eff(
  model = corr.brm,
  raw_data = dat$corr,
  title = "Correlated interaction strengths") 

# check posterior predictive ability of model
pp_check(corr.brm, type="boxplot")
# this model also looks good

# print out model summary
# (see MS_functions.R for more details)
corr.summary <- model_summary(corr.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/Correlated_model_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
    plot_prior_vs_posterior(corr.brm)
dev.off()

# remove model for RAM space
rm(corr.brm)

# scaled + correlated interaction strengths ####
s.c.brm <- readRDS( "data/Bayesian_stability_s_c.RDS")
# marginal efffects plot ####
s.c.marg.eff <- plot_marg_eff(
  model = s.c.brm,
  raw_data = dat$s.c,
  title = "Scaled & correlated") 

# check posterior predictive ability of model
pp_check(s.c.brm, type="boxplot")
# this model also looks good

# print out model summary
# (see MS_functions.R for more details)
s.c.summary <- model_summary(s.c.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/S_c_model_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
    plot_prior_vs_posterior(s.c.brm)
dev.off()

# remove model for RAM space
rm(s.c.brm)

# make table of model summaries
model.summaries <- rbind(rand.summary[[1]],
                        scale.summary[[1]],
                        corr.summary[[1]],
                        s.c.summary[[1]])
column_1 <- row.names(model.summaries)
model.summaries <- as.data.frame(model.summaries,
                                 row.names = NA)
model.summaries$variable <- column_1
model.summaries$model <- rep(c("Random",
                               "Scaled",
                               "Correlated",
                               "S.C"),
                             each = 5)
model.summaries <- model.summaries[,c(5,4,1:3)]

# write table as csv for MS
write.csv(model.summaries,
          "data/figures/stability_model_CrI.csv",
          row.names = FALSE)

# combine figures for MS
png(filename = "data/figures/stability_figure.png",
    width = 6.5,
    height = 9,
    units = "in",
    res = 300)
grid.arrange(
  rand.marg.eff +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    annotate(
      geom = "text",
      x = 5,
      y = 0.23,
      label = "A",
      size = 6
    ), 
  scale.marg.eff +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    annotate(geom = "text",
             x = 5,
             y = 0.09,
             label = "B",
             size = 6),
  corr.marg.eff +
    theme(legend.position = "none")+
    annotate(geom = "text",
             x = 5,
             y = 0.19,
             label = "C",
             size = 6),
  s.c.marg.eff +
    theme(legend.position = "none")+
    annotate(geom = "text",
             x = 5,
             y = 0.13,
             label = "D",
             size = 6),
  ncol = 2)
dev.off()


# stability ~ fw measures plots
dat <- dat.df %>% filter(.id == "random")
dat <- dat %>%
  mutate(num_0 = S**2 - L,
         prop_0 = num_0 / S**2)

ggplot(aes(x = log(num_0), y = log(stab), color = pca1), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")
ggplot(aes(x = prop_0, y = stab, color = pca1), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")


ggplot(aes(x = C, y = stab, color = pca1), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")

summary(lm(stab~prop_0, data = dat))

ggplot(aes(x = pca1, y = stab, color = pca1), data = dat)+
  geom_point(alpha = 0.2)


ggplot(aes(x = S, y = prop_0, color = pca1), data = dat)+
  geom_point(alpha = 0.2)

ggplot(aes(x = num_0, y = stab, color = pca1), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")

ggplot(aes(x = pca1, y = prop_0, color = pca1), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")

ggplot(aes(x = pca1, y = prop_0, color = stab), data = dat)+
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")
