# interpret fw measures bayesian model output
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(brms)
source("functions/MS_functions.R")

theme_set(theme_bw())

# simulation output data
# needed for points on marginal effects plot
dat <- readRDS("data/stability_results.RDS")
dat <- dat %>% 
  # select one model output e.g. random
  filter(.id == "random") %>%
  # calculate normalized G and V
  mutate(G = (L / (S* T+ S * I)) / S,
         V = (L / (S * B + S * I)) / S) %>%
  # only keeping some variables for Ram space
  select(Site, pca1, S, C, L, G, V, mean.TL, max.TL)

# Species ####
S.brm <- readRDS("data/Bayesian_FW_S.RDS")
S.marg.eff <- plot_marg_eff_fw(model = S.brm,
                               raw_data = dat,
                               raw_y = "S",
                               title = "No. Species") 

# check posterior predictive ability of model
pp_check(S.brm, type="boxplot")
# Boxplots largely overlap

# print out model summary
# (see MS_functions.R for more details)
S.summary <- model_summary(S.brm)

# Links ####
L.brm <- readRDS("data/Bayesian_FW_L.RDS")

# marginal efffects plot ####
L.marg.eff <- plot_marg_eff_fw(model = L.brm,
                            raw_data = dat,
                            raw_y = "L",
                            title = "No. Links") 

# check posterior predictive ability of model
pp_check(L.brm, type="boxplot")
# Boxplots largely overlap

# print out model summary
# (see MS_functions.R for more details)
L.summary <- model_summary(L.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/L_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
plot_prior_vs_posterior(L.brm)
dev.off()

# remove model for RAM space
rm(L.brm)

#================================
# connectance ####
C.brm <- readRDS("data/Bayesian_FW_C.RDS")

# marginal efffects plot ####
C.marg.eff <- plot_marg_eff_fw(model = C.brm,
                               raw_data = dat,
                               raw_y = "C",
                               title = "Connectance") 

# check posterior predictive ability of model
pp_check(C.brm, type="boxplot")
# Predictive ability looks good

# print out model summary
# (see MS_functions.R for more details)
C.summary <- model_summary(C.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/C_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
plot_prior_vs_posterior(C.brm)
dev.off()

# remove model for RAM space
rm(C.brm)

#================================
# Generality ####
Gen.brm <- readRDS("data/Bayesian_FW_Gen.RDS")

# marginal efffects plot ####
Gen.marg.eff <- plot_marg_eff_fw(model = Gen.brm,
                               raw_data = dat,
                               raw_y = "G",
                               title = "Normalized Generality") 

# check posterior predictive ability of model
pp_check(Gen.brm, type="boxplot")
# looks good

# print out model summary
# (see MS_functions.R for more details)
Gen.summary <- model_summary(Gen.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/Gen_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
plot_prior_vs_posterior(Gen.brm)
dev.off()

# remove randon model for RAM space
rm(Gen.brm)

#================================
# Vulnerability ####
Vul.brm <- readRDS("data/Bayesian_FW_Vul.RDS")

# marginal efffects plot ####
Vul.marg.eff <- plot_marg_eff_fw(model = Vul.brm,
                                 raw_data = dat,
                                 raw_y = "V",
                                 title = "Normalized Vulnerability") 

# check posterior predictive ability of model
pp_check(Vul.brm, type="boxplot")
# looks good

# print out model summary
# (see MS_functions.R for more details)
Vul.summary <- model_summary(Vul.brm)

# plot prior vs posterior
# for supplemental
png(filename = "data/figures/supplemental/Vul_prior_vs_posterior.png",
    width = 190, height = 190, units = "mm", res =300)
plot_prior_vs_posterior(Vul.brm)
dev.off()

# remove randon model for RAM space
rm(Vul.brm)


# make table of model summaries
model.summaries <- rbind(L.summary[[1]],
                         C.summary[[1]],
                         Gen.summary[[1]],
                         Vul.summary[[1]])
column_1 <- row.names(model.summaries)
model.summaries <- as.data.frame(model.summaries,
                                 row.names = NA)
model.summaries$variable <- column_1
model.summaries$model <- rep(c("L",
                               "C",
                               "SD Gen",
                               "SD Vul"),
                             each = 5)
model.summaries <- model.summaries[,c(5,4,1:3)]

# write table as csv for MS
write.csv(model.summaries,
          "data/figures/FW_model_CrI.csv",
          row.names = FALSE)

# combine figures for MS
png(filename = "data/figures/FW_figure.png",
    width = 6.5, height = 9, units = "in", res =300)
grid.arrange(L.marg.eff +
               theme(legend.position = "none",
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank())+
               annotate(geom = "text",
                        x = -3,
                        y = 175,
                        label = "A",
                        size = 6),
             C.marg.eff +
               theme(legend.position = "none",
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank())+
               annotate(geom = "text",
                        x = -3,
                        y = 0.25,
                        label = "B",
                        size = 6),
             Gen.marg.eff +
               theme(legend.position = "none",
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank())+
                       annotate(geom = "text",
                                x = -3,
                                y = 0.25,
                                label = "C",
                                size = 6),
             Vul.marg.eff+
               theme(legend.position = "none")+
               annotate(geom = "text",
                        x = -3,
                        y = 0.4,
                        label = "D",
                        size = 6),
             ncol = 1)
dev.off()
