# bayesian fw measures

library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms)

theme_set(theme_bw())

#read in simulation data
dat.df <- readRDS("data/stability_results.RDS")
# only keeping some variables for Ram space
# and only from one interaction strength estimate e.g. "random"
dat.df <- dat.df %>% 
  filter(.id == "random") %>%
  select(Site, pca1, C, L, Gensd, Vulsd)

# read in model results
C.brm <- readRDS("data/Bayesian_FW_C.RDS")
L.brm <- readRDS("data/Bayesian_FW_L.RDS")
Gen.brm <- readRDS("data/Bayesian_FW_Gen.RDS")
Vul.brm <- readRDS( "data/Bayesian_FW_Vul.RDS")

plot_marg_eff <- function(model, raw_y, raw_data, title = NULL){
  marg_eff <- marginal_effects(model, method = "fitted")
  marg_eff_df <- as.data.frame(marg_eff$pca1)
  #return(marg_eff_df)
  ggplot()+
    geom_point(data = raw_data,
               aes(x = pca1,
                   y = raw_data[,raw_y],
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
    ylab(raw_y) +
    xlab("PCA1 (AMD Stress)") +
    ggtitle(title)
}

plot_marg_eff(model = L.brm,
              raw_data = dat.df,
              raw_y = "L",
              title = "Number of links") 
plot_marg_eff(model = C.brm,
              raw_data = dat.df,
              raw_y = "C",
              title = "Connectance") 
plot_marg_eff(model = Gen.brm,
              raw_data = dat.df,
              raw_y = "Gensd",
              title = "SD normalized Generality")
plot_marg_eff(model = Vul.brm,
              raw_data = dat.df,
              raw_y = "Vulsd",
              title = "SD normalized Vulnerability")
coef_prob_not_0 <- function(model, coeficient){
  post_samples <- posterior_samples(model)[,coeficient]
  probability <- sum(post_samples!=0) / length(post_samples)
  return(probability)
}
coef_prob_not_0(L.brm, "b_pca1")
coef_prob_not_0(L.brm, "b_Intercept")
coef_prob_not_0(C.brm, "b_pca1")
coef_prob_not_0(C.brm, "b_Intercept")
coef_prob_not_0(Gen.brm, "b_pca1")
coef_prob_not_0(Gen.brm, "b_Intercept")
coef_prob_not_0(Vul.brm, "b_pca1")
coef_prob_not_0(Vul.brm, "b_Intercept")


plot_posterior <- function (model, raw_data, raw_y, title = NULL, rel_min_height = 0){
  newdata <- raw_data %>% distinct(pca1, Site)
  fitted_new <- fitted(model,
                       newdata = newdata,
                       re_formula = NA,
                       summary = FALSE)
  colnames(fitted_new) <-  newdata %>% 
    pull(Site) #rename columns as the sites that correspond to pca1
  fitted_plot <- as.data.frame(fitted_new)
  
  fitted_plot %>% 
    gather(Site, raw_y) %>%
    merge(newdata) %>% 
    ggplot(aes(x = .[,"raw_y"],
               y = reorder(Site, pca1),
               fill = reorder(Site, pca1)))+
    geom_density_ridges(rel_min_height = rel_min_height)+
    ggtitle(paste("Posterior distribution of", raw_y, "as predicted from\n", title, "model"))+
    scale_fill_manual(values = colorRampPalette(rev(brewer.pal(
      n = 11,
      name = "Spectral")))(25)) +
    theme(legend.position = "NULL") +
    ylab("Site") +
    xlab(raw_y)
}

plot_posterior(L.brm, dat.df, "L", title = "L.brm")
plot_posterior(C.brm, dat.df, "C", title = "C.brm")
plot_posterior(Gen.brm, dat.df, "Gensd", title = "Gen.brm")
plot_posterior(Vul.brm, dat.df, "Vulsd", title = "Vul.brm")
