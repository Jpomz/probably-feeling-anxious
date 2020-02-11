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
#=======================================#
# Functions for figures and results ####
#=======================================#
# Marginal effects plot
plot_marg_eff <- function(model, raw_data, title = NULL){
  # marginal effects
  marg_eff <- marginal_effects(model, method = "fitted")
  # store as data frame
  marg_eff_df <- as.data.frame(marg_eff$pca1)
  
  # plot raw data as points
  # plot marginal effects line with 95% Credible Interval
  ggplot()+
    geom_point(data = raw_data,
               aes(x = pca1,
                   y = stab,
                   color = reorder(Site, -pca1)
               ),
               alpha = 0.3,
               position = position_jitter(width=0.1))+
    geom_ribbon(data = marg_eff_df,
                aes(x = pca1,
                    ymin = lower__,
                    ymax = upper__),
                alpha = 0.15)+
    geom_line(data = marg_eff_df,
              aes(x = pca1,
                  y = estimate__),
              color = "blue",
              size=1) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),
                                 title="Site"))+
    scale_color_manual(values = colorRampPalette((brewer.pal(
      n = 11,
      name = "Spectral")))(25))+
    ylab(expression(paste("Stability, ", italic("s")))) +
    xlab("PCA1") +
    ggtitle(title)
}

# slightly modified function for FW_measures
plot_marg_eff_fw <- function(model, raw_y, raw_data, title = NULL){
  marg_eff <- marginal_effects(model, method = "fitted")
  marg_eff_df <- as.data.frame(marg_eff$pca1)
  ggplot()+
    geom_point(data = raw_data,
               aes(x = pca1,
                   y = raw_data[,raw_y],
                   color = reorder(Site, -pca1)
               ),
               alpha = 0.3,
               position = position_jitter(width=0.1))+
    geom_ribbon(data = marg_eff_df,
                aes(x = pca1,
                    ymin = lower__,
                    ymax = upper__),
                alpha = 0.15)+
    geom_line(data = marg_eff_df,
              aes(x = pca1,
                  y = estimate__),
              color = "blue",
              size=1) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),
                                 title="Site"))+
    scale_color_manual(values = colorRampPalette((brewer.pal(
      n = 11,
      name = "Spectral")))(25))+
    ylab(raw_y) +
    xlab("PCA1") +
    ggtitle(title)
}

# plot prior and posterior distributions
# this function assumes priors are the same for 
# both parameters
plot_prior_vs_posterior <- function(model,
                                    prior_mu = 0,
                                    prior_sigma = 1){
  int_coef <- fixef(model)[1,1:2]
  slope_coef <- fixef(model)[2,1:2]
  
  plot.dat <- data.frame(
    value = c(rnorm(1000, prior_mu, prior_sigma), 
              rnorm(1000, int_coef[1], int_coef[2]),
              rnorm(1000, slope_coef[1], slope_coef[2])),
    sample = rep(c("Both Priors",
                   "Intercept posterior",
                   "Slope posterior"),
                 each = 1000))
  ggplot(plot.dat, aes(x = value, fill = sample))+
    geom_density(alpha = 0.75)
}
# summarize bayesian model results
model_summary <- function(model){
  
  # extract posterior samples
  post_model <- posterior_samples(model)
  
  # coefficient information
  int_coef <- quantile(exp(post_model$b_Intercept),
                       probs=c(0.025,0.5,0.975))
  slope_coef <- quantile(exp(post_model$b_pca1),
                         probs=c(0.025,0.5,0.975))
  prob_slope_declines <- sum(exp(post_model$b_pca1) < 1) /
    nrow(post_model)
  
  # Relative change across gradient calculations
  reference = -3.5599
  min_impacts = -0.8089
  max_impacts = 6.1127
  
  # calculate fitted response variable at 3 points
  ref_stab <- exp(post_model$b_Intercept +
                    post_model$b_pca1 * 
                    reference)
  min_impacts_stab <- exp(post_model$b_Intercept +
                            post_model$b_pca1 * 
                            min_impacts)
  max_impacts_stab <- exp(post_model$b_Intercept +
                            post_model$b_pca1 * 
                            max_impacts)
  
  # fold change with increasing impacts
  # ref to min_impacts
  ref_to_min <- quantile(ref_stab/min_impacts_stab,
                         probs = c(0.025, 0.5, 0.975))
  # ref to max impacts
  ref_to_max <- quantile(ref_stab/max_impacts_stab,
                         probs = c(0.025, 0.5, 0.975))
  # min to max impacts
  min_to_max <- quantile(min_impacts_stab/max_impacts_stab,
                         probs = c(0.025, 0.5, 0.975))
  datout <-rbind(int_coef,
                 slope_coef,
                 ref_to_min,
                 ref_to_max,
                 min_to_max)
  datout <- list(CrI = datout,
                 Prob_declines = prob_slope_declines)
  return(datout)
}
