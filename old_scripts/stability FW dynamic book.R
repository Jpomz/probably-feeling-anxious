# interaction strength 
# FW dynamic book chapter 4
# J Pomz
# April 2018

# calculate coefficients of Jacobian matrix

# cij = aij * Bistar
# aij = -(Fij / Bi*Bj)
# Fij = fj * Bj * pij
# fj = 0.0096 Wi **0.82
# pij = Xi / sum(X..) where X.. = all prey of j

# j = predator
# i = prey
# so aij = effects of predator j on prey i
# Fij = per capita ingestion rate of j on i
# fj = individual ingestion rate

library(plyr)
library(tidyverse)
library(Matrix)

# bernoulli trials
dat <- readRDS("data/500 bernouli trials.rds")
# metadata of AMD sites
amd <- readRDS("data/AMD fish invert dw abundance.RDS")
# split into list to match structure of dat
amd.list <- split(amd, list(amd$site))
# arrange by increasing dw
amd.list <- map(amd.list, arrange, avg.dw)
# list of vectors of dw
dw.list <- map(amd.list, pull, "avg.dw")
dw.list <- map(dw.list, ~.x*1000)
# equilibrium density ####
#estimate equilibrium density
xistar.list <- map(dw.list, ~3*(.x)**-0.98)
# equilibrium biomass ####
# list of vectors of density
# density is definitely calculated wrong...
#density.list <- map(amd.list, pull, "density")
Bi.list <- map2(xistar.list, dw.list, ~.x*.y) 

# # pull out one adjacency matrix to figure out functions
# M <- dat[[13]][[1]]
# dw <- dw.list[[13]]
#density <- density.list[[13]]

# get_pij() ####
# estimate predator preference of prey
get_pij <- function(A, x){
  #preds <- which(colSums(A)>0, arr.ind = TRUE)
  preys <- which(A ==1, arr.ind = TRUE)
  n.prey <- colSums(A)[unique(preys[,2])]
  x.prey <- x[preys[,1]]
  tot.prey.x <- NULL 
  for (i in unique(preys[,2])){
    index <- preys[which(preys[,2] == i),]
    if(class(index) == "integer"){
      index <- index[1]
    }
    if(class(index) == "matrix"){
      index <- index[,1]
    }
    tot.prey.x[[i]] <- sum(x[index])
  }
  prey.proportion <- x.prey /
      rep(tot.prey.x[unique(preys[,2])], n.prey)
  pij <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for(i in 1:length(prey.proportion)){
    pij[preys[i,][1], preys[i,][2]] <- prey.proportion[i]
  }
  return(pij)
}

# pull out first bernoulli trial from each site to figure out functions
test <- lapply(dat, "[[", 10)

# pij ####
pij.list <- map2(test, xistar.list,
                 possibly(get_pij, otherwise = NA))

# fj ####
fj.list <- map(dw.list, ~0.0096*.x**0.82)

# Fij ###
# Fij = fj * Bj, * pij
get_Fij <- function(B, f, pij){
  # B = vector of equilibrium biomass
  # f = predator feeding rate
  # pij = matrix of predator (col) preference of prey (row)
  index <- which(pij !=0, arr.ind = TRUE)
  Fij <- matrix(0, nrow = nrow(pij), ncol = ncol(pij))
  for(r in 1:nrow(index)){
    pred = index[r,][2]
    prey = index[r,][1]
    Fij[prey, pred] <- f[pred]*B[pred]*pij[prey, pred]
  }
  return(Fij)
}

Fij.list <- pmap(list(B = Bi.list,
                      f = fj.list,
                      pij = pij.list),
                 possibly(get_Fij, otherwise = NA))

# aij ####
get_aij <- function(Fij, B){
  index <- which(Fij !=0, arr.ind = TRUE)
  Aij <- matrix(0, nrow = nrow(Fij), ncol = ncol(Fij))
  for (r in 1:nrow(index)){
    pred = index[r,][2]
    prey = index[r,][1]
    Aij[prey, pred] <- -Fij[prey, pred] / 
      (B[pred]*B[prey])
  }
  # # set diagonal values
  # # basal species = -1
  # # non-basal = 1e-06
  # for(c in 1:ncol(Aij)){
  #   if (sum(Aij[,c])!=0){
  #     Aij[c,c] <- 1e-06
  #   }
  #   if (sum(Aij[,c])==0){
  #     Aij[c,c] <- -1
  #   }
  # }
  diag(Aij) <- 0
  return(Aij)
}

Aij.list <- map2(Fij.list, Bi.list, 
                 possibly(get_aij, otherwise = NA))


# diagonal matrix B
B <- llply(Bi.list, function (x){
  mat = matrix(0, nrow = length(x),
               ncol = length(x))
  diag(mat) <- x
  return(mat)
})

# C = B %*% Aij
C <- map2(B, Aij.list, 
          possibly(~.x%*%.y, otherwise = NA))

# eigen
e.vals <- map(C, possibly(eigen, otherwise = NA))
#e.vals[[1]][[1]]

# prop.neg <- llply(e.vals, function (x){
#   sum(Re(x[[1]])<0) / length(x[[1]])
# })


dom.eig <- llply(e.vals, function (x){
  Re(x[[1]][1])
})

llply(e.vals, function (x){
  x[[1]][1]
})

dom.eig.df <- ldply(dom.eig)
pca.axis <- readRDS("data/pca axis 26 sites.rds")[,c(4, 2)]
dom.eig.df <- left_join(dom.eig.df, pca.axis, by = c(".id" = "site"))

ggplot(dom.eig.df, aes(x = pca1, y = V1)) +
  geom_point() +
  stat_smooth(method = "lm")# +
#  ylim(-3.5, -2.5)

prop.neg.df <- ldply(prop.neg)
prop.neg.df <- left_join(prop.neg.df, pca.axis, by = c(".id" = "site"))
ggplot(prop.neg.df, aes(x = pca1, y = V1)) +
  geom_point()



# from Sauve et al 2016 plant paper
stability <- function(J, s2){
  # J is the community matrix  (i.e. the jacobian matrix, see functions jacobian_binary() or jacobian_quant()) with a null diagonal, dim(J) = S*S
  # s2 is an arbitrary parameter such that the real part of the greatest eigenvalue of J-s2*I is negative
  
  test1 <- (dim(J)[1] == dim(J)[2]) # Is J a square matrix?
  test2 <- FALSE
  if (test1 == TRUE){
    S <- dim(J)[1]
    test2 <- which(diag(J) != vector("numeric", S)) # Does J have a null diagonal?
  }
  
  if ((test1 == TRUE) & (length(test2) == 0)){ # if J is a square matrix with a null diagonal
    S <- dim(J)[1] # S is the number of species in the network
    s1 <- 0
    I <- diag(S)
     E1 <- max(Re(eigen(J-s1*I, only.values = T)$values))
     E2 <- max(Re(eigen(J-s2*I, only.values = T)$values))
    # E1 <- (Re(eigen(J-s1*I, only.values = T)$values))[1]
    # E2 <- (Re(eigen(J-s2*I, only.values = T)$values))[1]
    if ((E1>=0) & (E2<0)){ # if s2 is well chosen and the system is not already stable
      while ((s2-s1)>=10^-4){
        stab <-(s1+s2)/2
        E1 <-max(Re(eigen(J-stab*I, only.values = T)$values))
        if (E1>=0){
          s1<-stab
        }
        else {
          s2<-stab
        }
      }
      return(stab)
    }
    
    if (E1<0){
      stop("J corresponds to a stable system.")
    }
    if (E2>=0){
      stop("s2 is not high enough.")
    }
  }
  else { # if J is not a square matrix with a null diagonal
    if (test1 == FALSE){
      stop("J is not a square matrix.")
    }
    if (length(test2) > 0){
      stop("J does not have a null diagonal.")
    }
  }
}

stab <- map(C, possibly(stability, otherwise = NA), s2 = 1)
stab.df <- ldply(stab)
stab.df <- left_join(stab.df, pca.axis, by = c(".id" = "site"))
ggplot(stab.df, aes(x = pca1, y = V1)) +
  geom_point() +
  stat_smooth(method = "lm")


test2 <- llply(test, function (x) {
  diag(x) = 0
  x[lower.tri(x)] = 0;x
})
set.seed(123)
J <- map(test2, jacobian_binary)
stab_sauve <- ldply(map(J, stability, s2 = 1))
stab_sauve <- left_join(stab_sauve, pca.axis, by = c(".id" = "site"))
ggplot(stab_sauve, aes(x = pca1, y = V1)) +
  geom_point() +
  stat_smooth(method = "lm")
