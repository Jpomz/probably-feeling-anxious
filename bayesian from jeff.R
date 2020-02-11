# Bayesian attempt
library(tidyverse)
library(RColorBrewer)
library(ggridges)
library(brms) 

theme_set(theme_bw())
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
dat$id <- dat$.id
ggplot(dat, aes(x = stab,
                   y = reorder(Site,pca1),
                   fill = reorder(Site,pca1))) +
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
  facet_wrap(~id)+
  coord_cartesian(xlim = c(0, 0.11))
# sites increas in AMD stress bottom to top
# stability increases R to L
# e.g. impacted sites "more" stable


# Bayesian model with pca1 as fixed effect and site as random effect [Is this correct syntax?]: YEP. Looks good.
rand.brm <- brm(stab ~ pca1 + (1|Site),
                 data = rand.j,
                 family = Gamma(link = "log"),
                 prior=c(prior(normal(0, 1),class = "b"),
                         prior(normal(0, 1),class = "Intercept")),
                 chains = 2,
                 iter = 1000,
                 cores = 4)

print(rand.brm,prior=TRUE)
# interpret output:
# what does the group-level effects mean? - Each site will have a different intercept, and those intercepts together will have a mean and standard deviiation from the overall intercept. This is the average "deviation' of each site's mean from the model intercept. I don't know why the word "family" is used. Not very helpful.

# pop level:
# pca1 = -0.27 --> with every 1 unit increase in pca1 gradient, the stability metric decreases by 0.27 [CrI -0.43, -0.14] ??
#----Not quite: The mean of y is multiplied by exp(-0.27) = 0.75 for each unit increase in x. This is why it looks non-linear
#----on the graph. e.g. if stability is 0.8 when x = -3, then stability would become 0.8*exp(-0.27)=0.8*0.75=0.6 when x=-2, and 
#----0.45 when x = -1, and so on. That allows y to never get below zero (or become equal to zero)

# family specific:
# no idea what this means -  IT'S THE SECOND PARAMETER OF THE GAMMA DISTRIBUTION. NOT REALLY IMPORTANT HERE.


plot(rand.brm)
# are the left panels the posterior distributions of the parameter? YES

marginal_effects(rand.brm, method="fitted")
# this looks to be a quadratic response
# is there any way to test that? or does this suggest we should actually look at a non-linear brm() e.g. stab ~ pca1 + I(pca1^2) + (1|Site)

# Gamma will always look curved b/c of the multiplicative nature of the log-link. However, it will look linear when the y-axis is plotted on the log scale. You can 
# try that with the ggplot code below. I've muted it for now, but if you make scale_y_log10(), then the line magicaly becomes linear (but won't go below zero, so it looks weird when numbers are near zero)




#------remake the marginal effects plot, and add in the raw data.
marg_eff <- marginal_effects(rand.brm, method = "fitted") #this gets the values to plot
plot_marg_eff <- as.data.frame(marg_eff$pca1) # data frame with predicted stability and CrI across pca1 gradient

ggplot()+
  geom_point(data = rand.j,
             aes(x = pca1,
                 y = stab,
                 color = reorder(Site, -pca1)),
             alpha = 0.2,
             position = position_jitter(width=0.1))+
  #scale_y_log10()+ #This will be linear on the log scale, but it will look weird b/c of lots of values near zero
  geom_ribbon(data = plot_marg_eff,
              aes(x = pca1,
                  ymin = lower__,
                  ymax = upper__),
              alpha = 0.25)+
  geom_line(data = plot_marg_eff,
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
  xlab("PCA1 (AMD Stress)")

# no Site level info? For this model, it has site levels in the random effects. We can simulate data from the posterior and plot the posterior
# distribution for each site using the code below.

newdata <- rand.j %>% distinct(pca1, Site) 

# For each unique value of pca1 (n=25) in newdata, return the posterior distribution of stability from the rand.brm model
fitted_new <- fitted(rand.brm,
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
  geom_density_ridges()+
  ggtitle("Posterior distribution of stability as predicted from rand.brm model")+
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(
    n = 11,
    name = "Spectral")))(25)) +
  theme(legend.position = "NULL") +
  ylab("Site") +
  xlab(expression(paste("Stability, ", italic("s"))))





# from here I get a bit lost and not sure where to go. 
# having joyplots of posterior dists would be cool (similar to figure at beginning of script), but not necessary if there is a better way. 


