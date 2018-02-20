# stability sensu Tang et al. 2012

# 1) remove lower triangle
# 2) Parameterize M
# 3) Calculate stability criteria (d == 0)

library(plyr)
library(tidyverse)

dat <- readRDS("data/100 bernouli trials 2018-02-14.rds")

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
dw <- mutate(dw, kg.m2 = (avg.dw / 1000) * 16.66667)
dw <- split(dw, list(dw$site))

# xistar ####
# could also calc as: xi = kg.m2 * density
xistar <- llply(dw, function (x){
  (x$avg.dw / 1000) * x$density
})

# # equilibrium biomass
# # xistar = 10^(x0 + 3 * gamma + epsilon) * mi^(1 + gamma)
# # mi = body mass in kg (avg.dw = g, so divide by 1000 to = kg)
# xistar <- llply(dw, function (x){
#   kg = x$avg.dw / 1000
#   x0 = -1.16
#   gamma = -0.675
#   epsilon = rnorm(n = 1, mean = 0, sd = 0.1)
#   result <- 10^(x0 + 3 * gamma + epsilon) * kg^(1 + gamma)
#   return(result)
# })

# aij ####
# search rate
# aij = 10^a0 * mi^bi * f(kij)
# a0 = -3.50
# mi = body mass (kg m^-2) = x$kg.m2
# f(kij) = function of kij
# kij = mj / mi
# f(kij) = kij^0.46 / (1 + kij^kappa)
# kappa = 2
# bi = N(-0.15, 0.0025); rnorm(1, -0.15, 0.0025)
# bi = constant ####
# setting bi to -0.15 for now!!!!!!!

# kij ####
# kij = mj / mi
# avg.dw = mg / 0.06 m^2
# (avg.dw / 1000) * 16.66667 ==> kg / m^2
# doesn't matter as long as use same for both mi and mj
# e.g. mg, kg*m^-2, etc, since proportion all units cancel anyways
kij <- map(dw, ~outer(X = .x$avg.dw, #.x$kg.m2,
                      Y = .x$avg.dw, #.x$kg.m2,
                      FUN = "/" ))
# f(kij) = kij^0.46 / (1 + kij^kappa)
f.kij <- map(kij, ~.x^0.46 / (1 + .x^2))

# # aij
# aij <- llply(kij, function (x){
#   matrix(0, ncol = ncol(x), nrow = nrow(x))
# })
# for(web in 1:length(dw)){
#   for(resource in 1:nrow(dw[[web]])){
#     for(consumer in 1:nrow(dw[[web]])){
#       #cons search rate for res
#       aij[[web]][resource, consumer] = 
#       10^-3.50 *
#       xistar[[web]][consumer]^-0.15 * 
#         # consumer body mass (kg*m^-2)
#       f.kij[[web]][resource,consumer] # f(kij)
#     }
#   }
# }

get_aij <- function (resource, consumer){
  aij = 10^-3.50 * consumer ^-0.15 *((resource / consumer)^0.46 / (1 + (resource / consumer)^2))
  return(aij)
}
xi.pairs <- map(xistar, ~expand.grid(.x, .x))
aij <- map(xi.pairs, ~get_aij(resource = .x[,1], consumer = .x[,2]))
aij <- map(aij, ~matrix(.x, nrow = sqrt(length(.x))))


dw.pairs <- map(dw, ~expand.grid(.x$avg.dw, .x$avg.dw))
aij.dw <- map(dw.pairs, ~get_aij(resource = .x[,1], consumer = .x[,2]))
aij.dw <- map(aij.dw, ~matrix(.x, nrow = sqrt(length(.x))))
# eij ####
# conversion efficiency, eij 
# herbivore eij = runif(n, min = 0.1, max = 0.3)
# carnivore eij = runif(n, min = 0.4, max = 0.6)
eij.fun <- function (n=1){
  # setting to 0.5 for now ####
  runif(n, min = 0.5, max = 0.5) #runif(n, min = 0.4, max = 0.6)
}

# interaction strengths
# mij = effect of prey resource on consumer
# mij = eij * aij * xistar
# mij[resource, consumer] = eij.fun * # conversion efficiency
#   aij[resource, consumer] * # search rate 
#   xistar[consumer]

# mji = effect of consumer on resource
# mji = -aij * xistar
# mij[consumer, resource] = 
#   -aij[resource, consumer] *
#   xistar[resource]

# aij
mij <- llply(kij, function (x){
  matrix(0, ncol = ncol(x), nrow = nrow(x))
})
for(web in 1:length(dw)){
  for(resource in 1:nrow(dw[[web]])){
    for(consumer in 1:nrow(dw[[web]])){
      if(resource < consumer){
        mij[[web]][resource, consumer] = 
          eij.fun() * # conversion efficiency
          aij[[web]][resource, consumer] * # search rate 
          xistar[[web]][consumer]
      }
      if (resource > consumer){
        mij[[web]][resource, consumer] = 
          -aij[[web]][resource, consumer] * # search rate 
          xistar[[web]][resource]
      }
    }
  }
}

# mij * bernouli trial
for(web in 1:length(dat)){
  for(trial in 1:length(dat[[web]])){
    dat[[web]][[trial]] <- mij[[web]]* 
      dat[[web]][[trial]]
  }
}

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
  out <- as.vector(out[out !=0])
  E2 = mean(out)
  stable = (sqrt(S*V)*(1+(E2 - E^2)/V)) - E < d
  rho = (E2 - E^2)/V
  result = list(stable = stable, rho = rho, vars = data.frame(d = d, S = S, E = E, V = V, E2 = E2))
  return(result)
}
stab_criterion(dat[[1]][[1]])

transpose(map(dat$Kiwi, stab_criterion)) %>%
  .$rho %>% flatten_dbl %>% sort

transpose(map(dat[[1]], stab_criterion)) %>%
  .$stable %>% flatten_lgl() %>% sum

pairs <- NULL # empty vector for pairwise interaction strengths
for(row in 1:nrow(M)){
  for(col in 1:ncol(M)){
    if (row < col){
      result = c(M[row, col], M[col, row])
      pairs <- rbind(result, pairs)
    }
  }
}

x <- map(dat[[25]], ~Re(eigen(.x)$values[1])) %>% flatten_dbl() 
plot(density(x))
density(x)


# stability ####
# start with stability (Tang analysis)
# sqrt(SV)(1+(E2 - E^2)/V) - E < d
# S = species numbers
# V = variance of off diagonal
# E = Mean of off diagonal
# E2 = mean product off diagonal pairs
# mean off diagonal
E = mean(mat[row(mat)!=col(mat)])
# variance off diagonal
V = sd(mat[row(mat)!=col(mat)])^2
# mean of off diagonal products
# mean(Mij*Mji)
up <- mat[upper.tri(mat)]
low <- mat[lower.tri(mat)]
E2 = mean(up * low)