# fw measure statistics
# S, L, C, B, I, T, Vulsd, Gensd
library(dplyr)
library(ggplot2)
library(robustbase)
library(lmtest)
library(forcats)

dat <- readRDS("data/stability_results.RDS") %>%
  filter(.id == "random")

# heteroskedasticity test format
# All response variables are heteroskedastic
bptest(S~pca1, data = dat)

# robust regressions
# testing for linear vs quadratic relationship
# anova result if p <0.05 then inclusion of additional terms (e.g. quadratic term) appropriate
anova(lmrob(S~pca1, data = dat),
      lmrob(S~pca1 +I(pca1^2), data = dat),
      test = "Wald")
S.mod <- lmrob(S~pca1 +I(pca1^2), data = dat)
summary(S.mod)

anova(lmrob(L~pca1, data = dat, 
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      lmrob(L~pca1 +I(pca1^2), data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      test = "Wald")
L.mod <- lmrob(L~pca1 +I(pca1^2),
               data = dat,
               control = lmrob.control(
                 maxit.scale = 700,
                 rel.tol = 1e-6))
summary(L.mod)

anova(lmrob(C~pca1, data = dat),
      lmrob(C~pca1 +I(pca1^2), data = dat),
      test = "Wald")
C.mod <- lmrob(C~pca1 +I(pca1^2), data = dat)
summary(C.mod)

anova(lmrob(B~1, data = dat),
      lmrob(B~pca1, data = dat),
      test = "Wald")
anova(lmrob(B~pca1, data = dat),
      lmrob(B~pca1 +I(pca1^2), data = dat),
      test = "Wald")
B.mod <- lmrob(B~pca1 +I(pca1^2), data = dat)
summary(B.mod)

anova(lmrob(I~pca1, data = dat),
      lmrob(I~pca1 +I(pca1^2), data = dat),
      test = "Wald")
I.mod <- lmrob(I~pca1 +I(pca1^2), data = dat)
summary(I.mod)

anova(lmrob(T~pca1, data = dat),
      lmrob(T~pca1 +I(pca1^2), data = dat),
      test = "Wald")
T.mod <- lmrob(T~pca1, data = dat)
summary(T.mod)

anova(lmrob(Vulsd~pca1, data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      lmrob(Vulsd~pca1 +I(pca1^2), data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      test = "Wald")
Vul.mod <- lmrob(Vulsd~pca1+I(pca1^2), data = dat,
                 control = lmrob.control(
                   maxit.scale = 700,
                   rel.tol = 1e-6))
summary(Vul.mod)

anova(lmrob(Gensd~pca1, data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      lmrob(Gensd~pca1 +I(pca1^2), data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-6)),
      test = "Wald")
Gen.mod <- lmrob(Gensd~pca1+I(pca1^2), data = dat,
               control = lmrob.control(
                 maxit.scale = 700,
                 rel.tol = 1e-6))
summary(Gen.mod)


# figures ####
fig_fun <- function(y,ylab,
                    quad = TRUE,
                    data = dat, x = "pca1",
                    maxit.scale = 200,
                    rel.tol = 1e-07,
                    linetype = "solid",
                    linecolor = "black"){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    geom_point(aes(color = fct_reorder(Site, pca1)))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_discrete(h.start = 80) +
    labs(x = "Mining gradient", y = ylab)+
    if(quad == TRUE)
      stat_smooth(method = "lmrob",
                  formula = y~x+I(x^2),
                  linetype = linetype,
                  color = linecolor, 
                  method.args = list(control =
                         lmrob.control(
                           maxit.scale = maxit.scale,
                           rel.tol = rel.tol))) 
    else
      stat_smooth(method = "lmrob",
                  linetype = linetype,
                  color = linecolor,
                  method.args = 
                    list(control =
                           lmrob.control(
                             maxit.scale = maxit.scale,
                             rel.tol = rel.tol))) 
}

fig_fun(y = "S", ylab = "No. Species")
ggsave("figures/S~gradient.png")
fig_fun(y = "L", ylab = "No. Links",
        maxit.scale = 700,
        rel.tol = 1e-6)
ggsave("figures/L~gradient.png")
fig_fun(y = "C", ylab = "Connectance")
ggsave("figures/C~gradient.png")
fig_fun(y = "B", ylab = "Prop. Basal",
        linetype = "dashed",
        linecolor = "red")
ggsave("figures/B~gradient.png")
fig_fun(y = "I", ylab = "Prop. Intermediate")
ggsave("figures/I~gradient.png")
fig_fun(y = "T", ylab = "Prop. Top", quad = FALSE)
ggsave("figures/T~gradient.png")
fig_fun(y = "Vulsd", ylab = "Vul SD",
        maxit.scale = 700, rel.tol = 1e-6)
ggsave("figures/Vul~gradient.png")
fig_fun(y = "Gensd", ylab = "Gen SD",
        maxit.scale = 700, rel.tol = 1e-6)
ggsave("figures/Gen~gradient.png")



