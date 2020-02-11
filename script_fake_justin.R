library(tidyverse)
library(ggridges)
library(brms) #Also need to install STAN from mc-stan.org. Not a trival task.





#-----------------------------------
#fake data for 200 runs on eight streams
#-----------------------------------
fake_stability_low <- data.frame(x = rgamma(800,.8,.5)/10, category = "low", id_run = rep(1:200,4),stream = as.factor(rep(1:4,each=200)))
fake_stability_high <- data.frame(x = rgamma(800,3.6,.5)/10, category = "high", id_run = rep(1:200,4),stream = as.factor(rep(5:8,each=200)))
fake_stability_all <- rbind(fake_stability_low,fake_stability_high)

# ggplot(fake_stability_all, aes(x=x,y=stream,fill=stream))+
#   geom_density_ridges(alpha=0.5)





#-----------------------------------
#bayesian Gamma GLMM of the fake data
#-----------------------------------
test_fake <- brm(x ~ stream + (1|id_run),data=fake_stability_all, family=Gamma(link="log"),
                 prior=c(prior(normal(0,1),class="b"),
                          prior(normal(0,1),class="Intercept")),
                 chains=1, iter=500,cores=4)

#model result
print(test_fake)


#-----------------------------------
#code to plot the model as a joy plot. Returns a plot of the posterior distribution of stability within each
#stream. It's bound at zero.
#-----------------------------------
marg_fake <-marginal_effects(test_fake,method="fitted")
marg_fake_fit <- fitted(test_fake, newdata=marg_fake$stream,summary=F)
columns_fake<- marg_fake$stream
colnames(marg_fake_fit) <- 1:8
marg_fake_fit <- as.data.frame(marg_fake_fit)
marg_fake_fit$iter <- 1:nrow(marg_fake_fit)

marg_fake_fit_plot <- as_tibble(marg_fake_fit) %>% 
  gather(stream, stability,-iter) %>% 
  ggplot(aes(x=stability/10,y=reorder(stream,-as.numeric(stream)),fill=stream))+
  geom_density_ridges()+
  coord_cartesian(xlim=c(0,0.1))+
  ylab("mining gradient")+
  xlab("stability")+
  ggtitle("Random Jij")+
  theme_classic()

marg_fake_fit_plot


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




