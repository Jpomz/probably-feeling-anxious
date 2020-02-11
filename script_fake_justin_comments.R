library(tidyverse)
library(ggridges)
library(brms) #Also need to install STAN from mc-stan.org. Not a trival task.
 

#-----------------------------------
#fake data for 200 runs on eight streams
#-----------------------------------
fake_stability_low <-
  data.frame(
    x = rgamma(800, .8, .5) / 10,
    category = "low",
    id_run = rep(1:200, 4),
    stream = as.factor(rep(1:4, each = 200))
  )
fake_stability_med <-
  data.frame(
    x = rgamma(800, 2.4, .5) / 10,
    category = "med",
    id_run = rep(1:200, 4),
    stream = as.factor(rep(5:8, each = 200))
  )
fake_stability_high <-
  data.frame(
    x = rgamma(800, 3.8, .5) / 10,
    category = "high",
    id_run = rep(1:200, 4),
    stream = as.factor(rep(9:12, each = 200))
  )
fake_stability_all <- rbind(fake_stability_low, fake_stability_med, fake_stability_high)

# imitating pca gradient values
fake_stability_all$gradient <- 
  scale(as.numeric(fake_stability_all$stream, scale = FALSE))

ggplot(fake_stability_all, aes(x=x,y=stream,fill=category))+
 geom_density_ridges(alpha=0.5)

#-----------------------------------
#bayesian Gamma GLMM of the fake data
#-----------------------------------

# just a continuous predictor
test_fake <-
  brm(
    x ~ gradient,
    data = fake_stability_all,
    family = Gamma(link = "log"),
    prior = c(prior(normal(0, 1), class = "b"),
              prior(normal(0, 1), class = "Intercept")),
    chains = 2,
    iter = 500,
    cores = 4
  )

#model result
print(test_fake)

#-----------------------------------
#code to plot the model as a joy plot. Returns a plot of the posterior distribution of stability within each
#stream. It's bound at zero.
#-----------------------------------
marg_fake <-marginal_effects(test_fake,method="fitted")
marg_fake_fit <- fitted(test_fake, newdata=marg_fake$gradient,summary=F)
columns_fake<- marg_fake$gradient
# colnames(marg_fake_fit) <- 1:8
marg_fake_fit <- as.data.frame(marg_fake_fit)
marg_fake_fit$iter <- 1:nrow(marg_fake_fit)
# 
# (marg_fake_fit_plot <- as_tibble(marg_fake_fit) %>% 
#   gather(stream, stability,-iter) %>% 
#   ggplot(aes(x=stability/10,y=reorder(stream,-as.numeric(stream)),fill=stream))+
#   geom_density_ridges()+
#   coord_cartesian(xlim=c(0,0.1))+
#   ylab("mining gradient")+
#   xlab("stability")+
#   ggtitle("Random Jij")+
#   theme_classic())

newdata <- marginal_effects(test_fake)$gradient
head(newdata)

ggplot(newdata, aes (y=estimate__, x = gradient)) +
  geom_point(data = fake_stability_all,
             aes(y = x, x = gradient),
             color = "gray",
             alpha = 0.5) +
  geom_line()+
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = "blue",
              alpha = 0.3) +
  scale_x_continuous("X") +
  scale_y_continuous("Y") +
  theme_bw()

#-----------------------------------
#summary stats of stability for each stream
#-----------------------------------

marg_fake_fit %>% 
  gather(stream, stability,-iter) %>% 
  group_by(stream) %>% 
  summarize(mean = mean(stability),
            sd = sd(stability),
            low95 = quantile(stability,probs=0.025),
            high95 = quantile(stability, probs=0.975)) %>% 
  mutate_if(is.numeric,round,2)


#interpretation of data - e.g. for stream 1 you could write: Given the model and data, I am 95% confident that #
#stability is between 0.13 and 0.18, with a mean value of 0.16 (Note from Jeff - your values will vary due to simulation error)


#-----------------------------------
# plots ####
#-----------------------------------
plot(marginal_effects(test_fake, points = TRUE))


# continuous + group level predictor
test_fake <-
  brm(
    x ~ gradient + (1|stream),
    data = fake_stability_all,
    family = Gamma(link = "log"),
    prior = c(prior(normal(0, 1), class = "b"),
              prior(normal(0, 1), class = "Intercept")),
    chains = 2,
    iter = 500,
    cores = 4
  )

#model result
print(test_fake)
plot(test_fake, pars = c("gradient", "stream"))
plot(marginal_effects(test_fake))$stream
marginal_effects(test_fake)$gradient
