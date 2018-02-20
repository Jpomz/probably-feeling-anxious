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

# aij ####
# function to calculate mass specific search rate
# should resource and consumer be kg/m2 (e.g. xistar), or individual bodymass (e.g. dw$avg.dw)?

get_aij <- function (resource, consumer,
                     a0 = -3.50, bi = -0.15){
  kij = resource / consumer
  aij = 10^a0 * consumer^bi *
    ((kij)^0.46 / (1 + (kij)^2))
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

# get all dw pairs
dw.pairs.kg.m2 <- map(dw, ~expand.grid(.x$kg.m2,
                                 .x$kg.m2))
# calculate aij for all dw pairs
aij.kg.m2 <- map(dw.pairs.kg.m2, 
              ~get_aij(resource = .x[,1],
                       consumer = .x[,2]))
# make aij.dw into matrices
aij.kg.m2 <- map(aij.kg.m2, ~matrix(.x,
                    nrow = sqrt(length(.x))))

# avg.dw / 1000
dw.pairs.1000 <- map(dw,
                    ~expand.grid(.x$avg.dw / 1000,
                                 .x$avg.dw / 1000))
aij.dw.1000 <- map(dw.pairs.1000, 
                 ~get_aij(resource = .x[,1],
                          consumer = .x[,2]))
aij.dw.1000 <- map(aij.dw.1000, ~matrix(.x,
                                    nrow = sqrt(length(.x))))
# eij ####
# conversion efficiency, eij 
# herbivore eij = runif(n, min = 0.1, max = 0.3)
# carnivore eij = runif(n, min = 0.4, max = 0.6)
eij.fun <- function (n=1){
  # setting to 0.5 for now ####
  runif(n, min = 0.5, max = 0.5) #runif(n, min = 0.4, max = 0.6)
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
          aij.kg.m2[[web]][resource, consumer] * 
          # search rate 
          xistar[[web]][consumer]
      }
      if (resource < consumer){
        mij[[web]][resource, consumer] = 
          -aij.kg.m2[[web]][resource, consumer] * 
          # search rate 
          xistar[[web]][resource]
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
  stable = (sqrt(S*V)*(1+(E2 - E^2)/V)) - E < d
  rho = (E2 - E^2)/V
  result = list(stable = stable, rho = rho, vars = data.frame(d = d, S = S, E = E, V = V, E2 = E2))
  return(result)
}

stab_criterion(M[[1]][[1]])

transpose(map(M$Kiwi, stab_criterion)) %>%
  .$rho %>% flatten_dbl %>% sort

transpose(map(M[[1]], stab_criterion)) %>%
  .$stable %>% flatten_lgl() %>% sum

 
x <- map(M[[1]], ~Re(eigen(.x)$values[1])) %>% flatten_dbl() 
plot(density(x))
abline(v = mean(x))

Re.eigen <- map(M, map, ~Re(eigen(.x)$values[1]))
Re.eigen <- lapply(Re.eigen, flatten_dbl)
plot(density(Re.eigen[[10]]))
