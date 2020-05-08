# estimating fish abundances with k-pass removal
# J Pomz

# this is based on sec 3.2 of the vignette: http://derekogle.com/fishR/examples/oldFishRVignettes/Depletion.pdf

# fishR is not available for R 3.6.0
# therefore, some functions were written to estimate the total number of fish based on three pass removal. 

# data with number of fish caught in three pass removal
fish <- read.csv("data/raw data/fish_k_pass_removal.csv")

# MLE function 
mle <- function(N0, X, k, T) { 
  (N0 + 0.5) * ((k * N0 - X - T)^k) - 
    (N0 - T + 0.5) * ((k * N0 - X)^k) 
}

# iteratively solve for N0 (total number of fish) based on 3 pass removal
N.fish <- function(vec){ 
  # vec = numeric vector where each value is number 
  # of fish caught in each successive pass
  T <-  sum(vec)
  k <-  length(vec)
  i <- seq(1, k)
  X <-  sum((k-i)*vec)
  T1 <- T
  mle.out <-  mle(N0 = T1, X = X, T = T, k = k)
  if(mle.out>0){
    while(mle.out > 0){
      T1 <- T1+1
      mle.out <- mle(N0 = T1, X = X, T = T, k = k)
    }
  }
return(T1 - T)
  #T1 = total estimated number of fish
  #T = total number of captured fish
  # T1 - T = number of fish to "add" to capture data
}


# 1st letter = site
# k = kiwi, b = Burke, i = italian, l = lankey, m = murray
(k.ang <- N.fish(vec = fish$Anguilla.australis[1:3])) # +2
(k.gob <- N.fish(vec = fish$Gobiomorphus.huttoni[1:3])) # +3
(k.fas <- N.fish(vec = fish$Galaxias.fasciatus[1:3])) # +0
(k.mac <- N.fish(vec = fish$Galaxias.maculatus[1:3])) # +2
(b.ang <- N.fish(vec = fish$Anguilla.australis[4:6])) # +0
(b.sal <- N.fish(vec = fish$Salmo.trutta[4:6]))# +0
(i.ang <- N.fish(vec = fish$Anguilla.australis[7:9]))# +0
(i.sal <- N.fish(vec = fish$Salmo.trutta[7:9]))# +0
(c.ang <- N.fish(vec = fish$Anguilla.australis[10:12]))# +0
(l.sal <- N.fish(vec = fish$Salmo.trutta[13:15]))# +0
(m.sal <- N.fish(vec = fish$Salmo.trutta[16:18]))# +0
