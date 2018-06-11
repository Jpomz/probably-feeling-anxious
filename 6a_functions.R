# 6a_functions
# jfpomeranz@gmail.com
# June 2018

# functions for for transforming probability matrices
# to adjacency matrices, and calculating interaction strengths

# b_trial() = bernoulli trial
# probability matrix to adjacency matrix
b_trial <- function (prob_matr){
  sum.out = 0
  # make sure that all matrices have at least one link
  while(sum.out == 0){
    out <- matrix(rbinom(length(prob_matr),
                         1,
                         prob = prob_matr),
                  ncol = ncol(prob_matr),
                  nrow = nrow(prob_matr))
    sum.out = sum(out[upper.tri(out)])
  }
  out
}
# function to set lower tri = 0
rm_low_tri <- function(adjacency){
  adjacency[lower.tri(adjacency, diag = TRUE)] <- 0
  adjacency
}
# function to estimate equilibrium biomass 
get_xistar <- function(N, M){
  (M / 1000) * N
}
# function to calculate search rate (aij)
# sensu tang et al. 2014
# add variability here?
get_aij <- function (xi, a0 = -3.50){
  xi.pairs = expand.grid(xi, xi)
  resource = xi.pairs[,1]
  consumer = xi.pairs[,2]
  bi = rnorm(1, -0.15, 0.05**2)
  kij = resource / consumer
  aij = 10^a0 * consumer^bi *
    (kij^0.46 / (1 + kij^2))
  aij = matrix(aij,
               nrow = length(xi),
               ncol = length(xi))
  return(aij)
}
# function to calculate conversion efficiency
# random number [0.4, 0.6], uniform distribution
eij <- function (n=1){
  runif(n, min = 0.4, max = 0.6)
}
# calculate mij (e.g. Jacobian matrix)
# each element is interaction strength
# upper triangle = negative (e.g. effect of pred on prey)
# lower triangle = positive (e.g. effect of prey on pred)
get_mij <- function(A, aij, xi){
  A.index = which(A > 0, arr.ind = TRUE)
  mij <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A.index)){
    mij[A.index[i,1], A.index[i,2]] = -aij[A.index[i,1], A.index[i,2]] *
      xi[A.index[i,1]]
    mij[A.index[i,2],A.index[i,1]] =eij() * 
      aij[A.index[i,1], A.index[i,2]] * xi[A.index[i,1]]
  }
  return(mij)
}