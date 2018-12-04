# stability sensu Tang et al. 2012

# 1) remove lower triangle
# 2) Parameterize M
# 3) Calculate stability criteria (d == 0)

library(plyr)
library(tidyverse)

dat <- readRDS("data/500 bernouli trials.rds")

# lower.tri ####
for(w in 1:length(dat)){
  for(t in 1:length(dat[[w]])){
    # make diagonal = 0
    diag(dat[[w]][[t]]) = 0 
    # make lower triangle = 0 
    dat[[w]][[t]][lower.tri(dat[[w]][[t]])] = 0
    # transpose upper to lower tri
    # this is for negative interaction strengths later
    dat[[w]][[t]] = dat[[w]][[t]] + t(dat[[w]][[t]])
  }
}


# parameterize M ####
set.seed(7968) # for reproducibility
dw <- readRDS("data/AMD fish invert dw abundance.RDS")
dw <- arrange(dw, site, avg.dw)
# avg.dw == grams, convert to kg 
dw <- mutate(dw, kg.m2 = (avg.dw / 1000) * density)
dw <- split(dw, list(dw$site))

# # xistar ####
xistar <- llply(dw, function (x){
  (x$avg.dw / 1000) * x$density
})

# # eastimate xistar as in Tang et al. 2014
# get_xistar <- function(dw){
#   x0 = -1.16
#   gamma = -.675
#   result = 10^(x0 + 3*gamma) * dw^(1+gamma)
#   return(result)
# }
# 
# xistar <- llply(dw, function (x){
#      get_xistar(x$avg.dw/1000)
#    })

# aij ####
# function to calculate mass specific search rate
# should resource and consumer be kg/m2 (e.g. xistar), or individual bodymass (e.g. dw$avg.dw)?

get_aij <- function (resource, consumer,
                     a0 = -3.50, bi = -0.15){
  kij = resource / consumer
  aij = 10^a0 * consumer^bi *
    (kij^0.46 / (1 + kij^2))
  return(aij)
}
# get all xistar pairs
xi.pairs <- map(xistar, ~expand.grid(.x, .x))
# calculate aij for all xistar pairs
aij <- map(xi.pairs, 
           ~get_aij(resource = .x[,1], 
                    consumer = .x[,2]))
# mak aij into matrices
aij <- map(aij, 
           ~matrix(.x, nrow = sqrt(length(.x))))


# eij ####
# conversion efficiency, eij 
# herbivore eij = runif(n, min = 0.1, max = 0.3)
# carnivore eij = runif(n, min = 0.4, max = 0.6)
eij.fun <- function (n=1){
  # setting to 0.5 for now ####
  runif(n, min = 0.5, max = 0.5) 
  #runif(n, min = 0.4, max = 0.6)
}

# interaction strengths
mij <- llply(aij, function (x){
  matrix(0, ncol = ncol(x), nrow = nrow(x))
})
for(web in 1:length(dw)){
  for(resource in 1:nrow(dw[[web]])){
    for(consumer in 1:nrow(dw[[web]])){
      if(resource > consumer){
        mij[[web]][resource, consumer] = 
          eij.fun() * # conversion efficiency
          aij[[web]][resource, consumer] * 
          # search rate 
          xistar[[web]][resource]
      }
      if (resource < consumer){
        mij[[web]][resource, consumer] = 
          -aij[[web]][resource, consumer] * 
          # search rate 
          xistar[[web]][consumer]
      }
    }
  }
}

# mij * bernouli trial
M <- llply(dat)

for(web in 1:length(dat)){
  for(trial in 1:length(dat[[web]])){
    M[[web]][[trial]] <- mij[[web]]* 
      dat[[web]][[trial]]
  }
}

 test <- lapply(M, "[[", 10)
# test <- llply(test, function (x){
#   x * 10e06
# })
# map(test, stability, s2 = 1)

# stability criterion ####
stab_criterion <- function(M){
  d = mean(diag(M))
  S = nrow(M)
  E = mean(M[row(M)!=col(M)])
  V = sd(M[row(M)!=col(M)])^2
  # mean of off diagonal products
  # mean(Mij*Mji)
  out <- NULL # empty vector for pairwise interaction strengths
  for(row in 1:nrow(M)){
    for(col in 1:ncol(M)){
      if (row < col){
        result = M[row, col] * M[col, row]
        out <- rbind(result, out)
      }
    }
  }
  out <- as.vector(out)
  E2 = mean(out)
  stable = sqrt(S*V)*(1+(E2 - E^2)/V) - E < d
  rho = (E2 - E^2)/V
  result = list(stable = stable, rho = rho, vars = data.frame(d = d, S = S, E = E, V = V, E2 = E2))
  return(result)
}

# real part leading eigen value 
re.eigen <- llply(M, function (x){
  map(x, ~Re(eigen(.x)$values[1])) %>%
    flatten_dbl()
})
re.eigen.df <- ldply(re.eigen) %>% 
  gather("trial", "re.eigen", -1)
pca.axis <- readRDS("data/pca axis 26 sites.rds")[,c(4, 2)]
re.eigen.df <- left_join(re.eigen.df, pca.axis, by = c(".id" = "site"))

# ggplot(re.eigen.df, aes(x = pca1, y = re.eigen))+
#   geom_point() +
#   geom_smooth(method = "lm") +
#   stat_summary(fun.y = "mean", color = "red",
#                geom = "point", size = 2, na.rm = TRUE) +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   coord_cartesian(ylim = (c(-2e-9, 5e-10))) +
#   theme_classic()
# 
# # figure out plot
# re.eigen.df %>%
#   filter(.id %in% c("Italia", "Hot")) %>%
# ggplot(aes(x = re.eigen))+
#   geom_density(adjust = 0.5) +
#   facet_wrap(~pca1, scales = "free")+
#   geom_vline(aes(xintercept = 0), linetype = 2) +
#   theme_classic() +
#   theme(strip.text = element_blank())
# 
# 
# ggplot(re.eigen.df, aes(x = re.eigen))+
#   geom_density(adjust = 2) +
#   facet_wrap(~pca1, scales = "free")+
#   geom_vline(aes(xintercept = 0), linetype =2,
#              size = 1) +
#   theme_classic() +
#   theme(strip.text = element_blank())

# proportion < 0
prop.stable <- re.eigen.df %>%
  mutate(stable = re.eigen < 0) %>%
  filter(stable == TRUE) %>%
  count(.id, pca1, stable) %>%
  mutate(prop = n / 500,
         sd.prop = sqrt((prop*(1 - prop))/500))

ggplot(prop.stable, aes(x = pca1, y = prop,
                        ymin = prop - sd.prop,
                        ymax = prop + sd.prop)) +
  geom_pointrange() +
  ylim(c(0, 1.0)) +
  #geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Mining gradient", y = "Proportion stable") +
  theme(axis.title = element_text(size = 20))
ggsave("figures/proportion stable.png")


prop.lm <- lm(prop~pca1, data = prop.stable)
summary(prop.lm)

normalized.re.eigen <- map(re.eigen,
          ~(.x /  mean(.x)))
normalized.re.eigen<- ldply(normalized.re.eigen) %>% 
  gather("trial", "re.eigen", -1)

normalized.re.eigen %>%
ggplot(aes(x = re.eigen, fill = .id))+
  geom_density(alpha = 0.5)+
  theme_classic() +
  ylim(c(0, 5)) +
  xlim(c(-1, 2))

# normalized.re.eigen %>%
#   filter(.id =="BV03") %>%
#   ggplot(aes(x = re.eigen, fill = .id)) +
#   geom_density() #+ xlim(c(-1, 1)) + ylim(c(0,10))

ldply(re.eigen) %>% 
  gather("trial", "re.eigen", -1) %>%
  ggplot(aes(x = re.eigen, fill = .id))+
  geom_density(alpha = 0.5#, adjust = 50
               )+
  theme_classic() +
  ylim(c(0, 1)) +
  geom_vline(xintercept = 0) #+
  #facet_wrap(~.id)

lapply(re.eigen, mean)


for(i in 1:length(re.eigen)){
  plot(density(re.eigen[[i]]))
  abline(v = 0)
}

llply(M, function (x){
  map(x, ~stab_criterion(.x))
})
 
# interaction correlations
get_rho <- function(M){
  E = mean(M[row(M)!=col(M)])
  V = sd(M[row(M)!=col(M)])^2
  # mean of off diagonal products
  # mean(Mij*Mji)
  out <- NULL # empty vector for pairwise interaction strengths
  for(row in 1:nrow(M)){
    for(col in 1:ncol(M)){
      if (row < col){
        result = M[row, col] * M[col, row]
        out <- rbind(result, out)
      }
    }
  }
  out <- as.vector(out)
  E2 = mean(out)
  rho = (E2 - E^2)/V
  return(rho)
}

rho <- llply(M, function (x){
  flatten_dbl(map(x, ~get_rho(.x)))
})
rho <- ldply(rho) %>%
  gather("trial", "rho", -1)
# this is from M-N manuscript (chapter 1)
# I need to re do this per reviewers comments from freshwater biology
# make sure to update with proper loadings in future. 
pca.axis <- readRDS("data/pca axis 26 sites.rds")[,c(4, 2)]
rho <- left_join(rho, pca.axis, by = c(".id" = "site"))

ggplot(rho, (aes(x = rho, fill = fct_reorder(.id, pca1))))+
  geom_density(alpha = 0.4) +
  theme_gray() +
  scale_fill_grey(start = 0, end = 1) +
  facet_wrap(~pca1) +
  theme(strip.text = element_blank())

rho %>%
  group_by(pca1) %>%
  summarize(mean = mean(rho), sd = sd(rho)) %>%
  ggplot(aes(x = pca1, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  stat_smooth(method = "lm")
rho.lm <- lm(rho~pca1, data = rho)
summary(rho.lm)
 
P <- M$Italia[[1]]
pair <- NULL # empty vector for pairwise interaction strengths
for(row in 1:nrow(P)){
 for(col in 1:ncol(P)){
   x = P[row, col] 
   y = P[col, row]
   out <- c(x, y)
   pair <- rbind(pair, out)
 }
}
plot(pair*1e6) 
abline(h = 0, v = 0)


test <- map(M$Italia, ~list(real = Re(eigen(.x)$values),
             im = Im((eigen(.x)$values)), 
             real.norm = Re(eigen(.x)$values) / 
               mean(Re(eigen(.x)$values)))) 

plot(test[[1]]$im ~ test[[1]]$real)
map(test, ldply)

get_pairs <- function (M){
  pair <- NULL
  for(row in 1:nrow(M)){
    for(col in 1:ncol(M)){
      if(M[row, col] !=0){
        x = M[row, col] 
        y = M[col, row]
        out <- c(x, y)
        pair <- rbind(pair, out)
      }
    }
  }
  pair
}

burke.pairs <- (map(M[[1]], get_pairs))
plot(distinct(burke.pairs) * 1e6)
abline(h = 0, v = 0)

italia.pairs <- map(M$Italia, get_pairs) 
plot(italia.pairs[[1]]*1e6)


x <- c()
y <- c()
for(row in 1:nrow(M)){
  for(col in 1:ncol(M)){
    if(M[row, col] !=0){
      x = M[row, col] 
      y = M[col, row]
      out <- c(x, y)
      pair <- rbind(pair, out)
    }
  }
}




# mean absolute value interaction strength
avg.intxn <- llply(M, function (x){
  flatten_dbl(map(x, ~mean(abs(.x))))
})

x <- ldply(map(avg.intxn, mean))
x <- left_join(x, pca.axis, by = c(".id" = "site"))
plot(V1~pca1, data = x)


ldply(avg.intxn) %>%
  ggplot(aes(x = V1, fill = .id))+
  geom_density(alpha = 0.5)




stab_value <- function(M){
  #d = mean(diag(M))
  S = nrow(M)
  E = mean(M[row(M)!=col(M)])
  V = sd(M[row(M)!=col(M)])^2
  # mean of off diagonal products
  # mean(Mij*Mji)
  out <- NULL # empty vector for pairwise interaction strengths
  for(row in 1:nrow(M)){
    for(col in 1:ncol(M)){
      if (row < col){
        result = M[row, col] * M[col, row]
        out <- rbind(result, out)
      }
    }
  }
  out <- as.vector(out)
  E2 = mean(out)
  stable = sqrt(S*V)*(1+(E2 - E^2)/V) - E 
  return(stable)
}

stab.val <- llply(M, function (x){
  flatten_dbl(map(x, stab_value))
})

ldply(stab.val) %>% gather("trial", "stab", -1) %>%
  group_by(.id) %>%
  summarize(stab = mean(stab)) %>%
  left_join(pca.axis, by = c(".id" = "site")) %>%
  ggplot(aes(x = pca1, y = stab)) +
  geom_point()


  ggplot(aes(x = stab, fill = .id))+
  geom_density(alpha = 0.5)
