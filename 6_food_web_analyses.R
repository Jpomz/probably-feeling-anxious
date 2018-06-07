# food web analyses

Pij <- readRDS("data/AMD_final_probability_matrices.RDS")


matr <- Pij[[1]]
set.seed(3049)

b_trial <- function (prob_matr){
  out <- matrix(rbinom(length(prob_matr),
                       1,
                       prob = prob_matr),
                ncol = ncol(prob_matr),
                nrow = nrow(prob_matr))
  out
}

set.seed(3049)
A <- b_trial(matr)

rm_low_tri <- function(adjacency){
  adjacency[lower.tri(adjacency, diag = TRUE)] <- 0
  adjacency
}

A <- rm_low_tri(adjacency = A)

taxa.dat <- readRDS("data/AMD_fish_invert_dw_abundance.RDS")
dat <- split(taxa.dat, list(taxa.dat$site))
dat <- dat[[1]]
dat <- dat[order(dat$avg.dw),]

get_xistar <- function(N, M){
  (M / 1000) * N
}

xistar <- get_xistar(N = dat$density, M = dat$avg.dw)

# add variability here?
get_aij <- function (xi, a0 = -3.50, bi = -0.15){
  xi.pairs = expand.grid(xi, xi)
  resource = xi.pairs[,1]
  consumer = xi.pairs[,2]
  kij = resource / consumer
  aij = 10^a0 * consumer^bi *
    (kij^0.46 / (1 + kij^2))
  aij = matrix(aij,
               nrow = length(xi),
               ncol = length(xi))
  return(aij)
}

aij <- get_aij(xistar)

eij <- function (n=1){
  runif(n, min = 0.4, max = 0.6)
}

get_mij <- function(A, aij){
  A.index = which(A > 0, arr.ind = TRUE)
  mij <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A.index)){
    mij[A.index[i,1],A.index[i,2]] = -aij[A.index[i,1], A.index[i,2]] *
      xi[A.index[i,1]]
    mij[A.index[i,2],A.index[i,1]] =eij() * 
      aij[A.index[i,1], A.index[i,2]] * xi[A.index[i,1]]
  }
  return(mij)
}

mij <- get_mij(A = A, aij = aij)

Re(eigen(mij)$values[1])
