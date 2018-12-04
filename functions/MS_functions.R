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


# function to convert to probability vector (vec) to probability matrix
get_prob_matr <- function(vec){
  m = matrix(vec, sqrt(length(vec)), sqrt(length(vec)))
}

# function to calculate relative abundance matrices
get_rel_ab <- function(vec, taxa){
  stopifnot(length(vec) == length(taxa))
  rel.ab <- vec / sum(vec)
  Nij <- matrix(0, length(vec), length(vec))
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  dimnames(Nij) <- list(taxa, taxa)
  Nij
}

# function to correct fish relative abundances
correct_fish_rel_ab <- function(matr, fish, cf = 10^3){
  if(any(fish %in% colnames(matr))){
    fish.cols <- which(colnames(matr) %in% fish)
    matr[,fish.cols] <- matr[,fish.cols] * cf
  }
  matr
}

# function to rescale variable x to [min, max]
scalexy <- function(x, min, max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

# function to plot heat map
# just using this for data exploration
plot_heat <- function(x, ...){
  heatmap.2(x,
            Rowv = NA,
            Colv = NA,
            scale = "none",
            trace = "none",
            dendrogram = "none",
            breaks = seq(0,1,0.01),
            key = F,
            labRow = NA,
            labCol = NA,
            main= "Probability of interaction",
            xlab = "Consumer",
            ylab = "Resource")
}

# function to set probabilities of forbidden taxa to 0
rm_niche <- function(matrix, taxa){
  for(name in (colnames(matrix)[colnames(matrix) %in% taxa])){
    matrix[,name] <- 0
  }
  matrix
}

# function to get fw_measures and dominant eigenvalue
get_measures <- function(matr, s2,
                         trials = 1,
                         scale.Jij = FALSE,
                         correlate.Jij = FALSE){
  result <- list()
  for(i in 1:trials){
    A <- rm_cycle(b_trial(matr))
    J <- jacobian_binary(A)
    if(scale.Jij == TRUE){
      J = scale.Jij(J)
    }
    if(correlate.Jij == TRUE){
      J = correlate.Jij(J)
    }
    stab = stability(J, s2 = s2)
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(stab = stab, fw_meas) 
  }
  return(result)
}

# function to scale interaction strengths by body size
scale.Jij <- function(J){
  # J is a Jacobian matrix where elements Jij = interaction strength
  # e.g. object from jacobian_binary()
  element.rank <- rank(1:nrow(J)) / nrow(J)
  scale.rank <- scalexy(element.rank, min = 0.25, max = 1.25)
  scale.grid <- expand.grid(rev(scale.rank), scale.rank)
  scale.vector <- scale.grid$Var1 * scale.grid$Var2
  scale.matrix <- matrix(scale.vector, nrow = nrow(J), ncol = ncol(J))
  #u.tri <- which(upper.tri(scale.matrix, diag = FALSE),
  #               arr.ind = TRUE)
  #l.tri <- cbind(u.tri[,2], u.tri[,1])
  #scale.matrix[l.tri] = scale.matrix[u.tri]
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

# function to construct food webs with given number of species (S) and connectance (C)
random_A <- function(S, C, ...){
  A <- 0
  while(sum(A) == 0){
    A <- matrix(rbinom(n = S^2, prob = C, size = 1), S, S)
    diag(A) <- 0
    rm_cycle(A)
  }
  return(A)
}

# function to get fw_measures and dominant eigenvalue for food webs with random structure
get_random_measures <- function(S,
                                C,
                                s2,
                                trials = 1,
                                scale.Jij = FALSE,
                                correlate.Jij = FALSE){
  result <- list()
  for(i in 1:trials){
    A <- random_A(S, C)
    J <- jacobian_binary(A)
    if(scale.Jij == TRUE){
      J = scale.Jij(J)
    }
    if(correlate.Jij == TRUE){
      J = correlate.Jij(J)
    }
    stab = stability(J, s2 = s2)
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(stab = stab, fw_meas) 
  }
  return(result)
}


# below this is trash?
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