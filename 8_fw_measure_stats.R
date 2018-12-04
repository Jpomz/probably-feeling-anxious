# fw measure statistics
# S, L, C, B, I, T, Vulsd, Gensd
library(dplyr)
library(ggplot2)
library(robustbase)
library(lmtest)
library(forcats)
library(gridExtra)

dat <- readRDS("data/stability_results.RDS") %>%
  filter(.id == "random")

# species stats
# filtering out one observation of each S:pca1 combination
S.dat <- dat %>%
  select(S, pca1, Site) %>%
  distinct
bptest(S~pca1, data = S.dat)
anova(lm(S~pca1, data = S.dat),
      lm(S~pca1 +I(pca1^2), data = S.dat))
AIC(lm(S~pca1, data = S.dat),
    lm(S~pca1 +I(pca1^2), data = S.dat))
# AIC gives both models equal weight
# going with 
S.mod <- lm(S~pca1 +I(pca1^2), data = S.dat)
summary(S.mod)

# heteroskedasticity test format
# All response variables are (Except Species no.) heteroskedastic
bptest(L~pca1, data = dat)

# robust regressions
# testing for linear vs quadratic relationship
# anova result if p <0.05 then inclusion of additional terms (e.g. quadratic term) appropriate
anova(lmrob(L~pca1, data = dat, 
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-5)),
      lmrob(L~pca1 +I(pca1^2), data = dat,
            control = lmrob.control(
              maxit.scale = 700,
              rel.tol = 1e-5)),
      test = "Wald")
L.mod <- lmrob(L~pca1 +I(pca1^2),
               data = dat,
               control = lmrob.control(
                 maxit.scale = 700,
                 rel.tol = 1e-5))
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
  ggplot(data, aes(x = -data[[x]], y = data[[y]])) +
    geom_point(aes(color = fct_reorder(Site, pca1)))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_discrete(h.start = 80) +
    labs(x = "AMD gradient", y = ylab)+
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

# fig 1A-C ####
ggplot(S.dat, aes(x = -pca1, y = S)) +
  geom_point(aes(color = fct_reorder(Site, pca1)))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_discrete(h.start = 80) +
  labs(x = "AMD gradient", y = "No. Species") +
  stat_smooth(method = "lm", 
              formula = y~x+I(x^2),
              color = "black") +
    annotate(geom = "text", x = -5.6, y = 33,
             label = "A", size = 6))

#ggsave("figures/S~gradient.png")
L.fig <- fig_fun(y = "L", ylab = "No. Links",
        maxit.scale = 700,
        rel.tol = 1e-5) +
  annotate(geom = "text", x = -5.6, y = 150,
           label = "B", size = 6)
#ggsave("figures/L~gradient.png")
C.fig <- fig_fun(y = "C", ylab = "Connectance") +
  annotate(geom = "text", x = -5.6, y = 0.3,
           label = "C", size = 6)
#ggsave("figures/C~gradient.png")

tiff(filename = "figures/fig_1.tiff",
     width = 190, height = 210, units = "mm", res =300)
grid.arrange(S.fig, L.fig, C.fig, ncol = 1)
dev.off()


fig_fun(y = "B", ylab = "Prop. Basal",
        linetype = "dashed",
        linecolor = "red")
#ggsave("figures/B~gradient.png")

# fig 2 A-D
I.fig <- fig_fun(y = "I", ylab = "Prop. Intermediate") +
  annotate(geom = "text", x = -5.6, y = 0.9,
           label = "A", size = 6)

T.fig <- fig_fun(y = "T", ylab = "Prop. Top", quad = FALSE)+
  annotate(geom = "text", x = -5.6, y = 0.75,
           label = "B", size = 6)

V.fig <- fig_fun(y = "Vulsd", ylab = "Vul SD",
        maxit.scale = 700, rel.tol = 1e-6)+
  annotate(geom = "text", x = -5.6, y = 0.17,
           label = "C", size = 6)

G.fig <- fig_fun(y = "Gensd", ylab = "Gen SD",
        maxit.scale = 700, rel.tol = 1e-6) +
  annotate(geom = "text", x = -5.6, y = 0.17,
           label = "D", size = 6)

tiff(filename = "figures/fig_2.tiff",
     width = 190, height = 190, units = "mm", res =300)
grid.arrange(I.fig, T.fig, V.fig, G.fig, ncol = 2)
dev.off()

