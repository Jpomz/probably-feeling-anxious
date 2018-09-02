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



scale.Jij <- function(J){
  # J is a Jacobian matrix where elements Jij = interaction strength
  # e.g. object from jacobian_binary()
  element.rank <- rank(1:nrow(J)) / nrow(J)
  scale.rank <- scalexy(element.rank, min = 0.25, max = 1.25)
  scale.grid <- expand.grid(rev(scale.rank), scale.rank)
  scale.vector <- scale.grid$Var1 * scale.grid$Var2
  scale.matrix <- matrix(scale.vector, nrow = nrow(J), ncol = ncol(J))
  u.tri <- which(upper.tri(scale.matrix, diag = FALSE),
                 arr.ind = TRUE)
  l.tri <- cbind(u.tri[,2], u.tri[,1])
  scale.matrix[l.tri] = scale.matrix[u.tri]
  scale.J <- scale.matrix * J
  return(scale.J)
}

scalexy <- function(x, min, max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

correlate.Jij <- function(J){
  # J is a JAcobian matrix where elements Jij = interaction strength
  # e.g. object from jacobian_binary()
  negative.index <- which(J < 0, arr.ind = TRUE)
  positive.index <- cbind(negative.index[,2], negative.index[,1])
  positive.strength <- abs(J[negative.index] * 0.7)
    #runif(n = nrow(negative.index), min = 0.6, max = 0.8))
  # 0.7 comes from Montoya et al. 2009
  J[positive.index] <- positive.strength
  return(J)
}



rm_cycle <- function(A, diag = TRUE){
  # if A eats B and B eats A, randomly remove one link
  # Default removes diagonal
  # to leave diagonal as is set diag = FALSE
  for(i in 1:nrow(A)){
  for(j in 1:ncol(A)){
    if(diag == FALSE){
      if(i == j){
        next
      }
    }
    if(A[i,j] == 1 & A[j,i] == 1){
      if(runif(1) < 0.5){
        A[i,j] <-  0
      }
      else{
        A[j,i] <- 0
      }
    }
  }
}
  return(A)
}


fwpointrange <- function(data, y, x = "pca1", ylab = NULL){
  ggplot(data, aes(x = data[[x]], y = data[[y]])) +
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x),
                 fun.ymax = function(x) mean(x) + sd(x),
                 geom = "pointrange") +
    theme_bw() +
    theme(axis.title = element_text(size = 20)) +
    if (is.null(ylab))
      labs(x = "Mining gradient", y = y) 
  else
    labs(x = "Mining gradient", y = ylab) 
}

joy.stability <- function(data, x = "stab",
                          xmin = -.09, xmax = 0.5,
                          title = NULL, scale = NULL,
                          bandwidth = NULL,
                          rel_min_height = NULL){
  require(ggridges)
  require(forcats)
  require(ggplot2)
  ggplot(data, aes(x = data[[x]],
                   y = fct_reorder(.id, pca1),
                   height = ..density.., 
                   fill = pca1)) +
    stat_density_ridges(scale = scale, rel_min_height = rel_min_height,
                        bandwidth = bandwidth, alpha = 0.8) +
    theme_ridges(center_axis_labels = TRUE) +
    labs(y = "mining gradient", x = "stability") +
    scale_fill_distiller(palette = "Spectral") +
    theme(legend.position = "NULL")+
    coord_cartesian(xlim = c(xmin, xmax)) +
    if (is.null(title))
      labs(title = NULL) 
  else
    labs(title = title)
}