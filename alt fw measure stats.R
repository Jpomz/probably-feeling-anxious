# trying gls() in nlme
library(nlme)
library(dplyr)
library(ggplot2)
library(mgcv)

dat <- readRDS("data/stability_results.RDS") %>%
  filter(.id == "random")

lm.mod <- lm(L~pca1+I(pca1^2)+I(pca1^3), data = dat)
dat$resi <- lm.mod$residuals
ggplot(data = dat, aes(y = resi, x = pca1))+
  geom_point()

gam_L <- gam(log10(L)~s(pca1), method = "REML", data = dat)
par(mfrow = c(2,2))
gam.check(gam_L)
summary(gam_L)

ggplot(dat, aes(y = L, x = pca1)) +
  geom_point()+
  stat_smooth(method = "gam", formula = y~s(x, bs = "ad", k = 9))
ggplot(dat, aes(y = C, x = pca1)) +
  geom_point()+
  stat_smooth(method = "gam", formula = y~s(x))
ggplot(dat, aes(y = Gensd, x = pca1)) +
  geom_point()+
  stat_smooth(method = "gam", formula = y~s(x))
ggplot(dat, aes(y = Vulsd, x = pca1)) +
  geom_point()+
  stat_smooth(method = "gam", formula = y~s(x))



l.mod <- gls(L~pca1, data = dat, weights = varPower())
l.mod2 <- gls(L~pca1+I(pca1^2), data = dat, weights = varPower())
l.mod3 <- gls(L~pca1+I(pca1^2)+I(pca1^3), data = dat, weights = varPower())
AIC(l.mod, l.mod2, l.mod3)

summary(l.mod)

plot(lm(L~pca1, data = dat))


power.mod <- gls(L~pca1, data = dat, weights = varPower())
power.mod2 <- gls(L~pca1+I(pca1^2), data = dat, weights = varPower())
power.mod3 <- gls(L~pca1+I(pca1^2)+I(pca1^3), data = dat, weights = varPower())
exp.mod <- gls(L~pca1, data = dat, weights = varExp())
fixed.mod <- gls(L~pca1, data = dat, weights = varFixed(~pca1))

plot(power.mod)
plot(power.mod2)
plot(power.mod3)

plot(gls(L~pca1, data = dat, weights = ~1/pca1))
plot(gls(L~pca1, data = dat, weights = ~pca1))
plot(gls(L~pca1+I(pca1^2), data = dat, weights = ~1/pca1))
plot(gls(L~pca1+I(pca1^2)+I(pca1^3), data = dat, weights = ~1/pca1))

plot(gls(log10(L)~pca1, data = dat, weights = NULL))
plot(gls(log10(L)~pca1, data = dat, weights = varPower()))
plot(gls(log10(L)~pca1, data = dat, weights = varExp()))
     

plot(gls(log10(L)~pca1+I(pca1^2), data = dat, weights = NULL))
plot(gls(log10(L)~pca1+I(pca1^2), data = dat, weights = varPower()))
plot(gls(log10(L)~pca1+I(pca1^2), data = dat, weights = varExp()))
