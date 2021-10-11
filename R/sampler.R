
# sampling from the p-value (calibration) rejection set union the knockoff rejection set
y_sampler_cond_Sj <- function(sample_size, calib_rej_reg, kn_rej_reg, j,
                              y.pack, X.pack){
  df <- y.pack$df
  y_Pi_Xnoj_res_norm2 <- y.pack$y_Pi_Xnoj_res_norm2[j]

  # mass of calibration rejection set
  rej_set <- calib_rej_reg
  if(!is.null(rej_set$left)){
    pval_rej_lb <- pt(rej_set$left, df = df, lower.tail = T)
    pval_rej_ub <- pt(rej_set$right, df = df, lower.tail = T)
    rej_mass <- sum(pval_rej_ub - pval_rej_lb)
  } else{
    rej_mass <- 0
  }

  # mass of (kn rejection set) - (calibration rejection set)
  rest_set <- interval_minus(kn_rej_reg, calib_rej_reg)
  if(!is.null(rest_set$left)){
    pval_rest_lb <- pt(rest_set$left, df = df, lower.tail = T)
    pval_rest_ub <- pt(rest_set$right, df = df, lower.tail = T)
    rest_mass <- sum(pval_rest_ub - pval_rest_lb)
  } else{
    rest_mass <- 0
  }

  total_mass <- rej_mass + rest_mass

  # can decide rejecting H_j or not based on the mass of the two sets
  if(total_mass == 0){
    return(0)
  }
  if(rest_mass == 0){
    return(1)
  }
  if(rej_mass == 0){
    return(-1)
  }

  sample_size <- 2 * ceiling(sample_size / 2)

  # sampling in the calibration rejection set and compute the
  # importance sampling weights
  rej_size <- sample_size / 2
  if(rej_size > 0){
    rej_p <- runif_intervals(pval_rej_lb, pval_rej_ub, rej_size)
    rej_t <- qt(rej_p, df = df, lower.tail = T)

    rej_weight <- rej_mass * 2
    rej_weights <- rep(rej_weight, rej_size)
  } else{
    rej_t <- NULL
    rej_weight <- 0
    rej_weights <- NULL
  }


  # sampling in (kn rejection set) - (calibration rejection set) and compute the
  # importance sampling weights
  rest_size <- sample_size - rej_size
  if(rest_size > 0){
    rest_p <- runif_intervals(pval_rest_lb, pval_rest_ub, rest_size)
    rest_t <- qt(rest_p, df = df, lower.tail = T)

    rest_weight <- rest_mass * 2
    rest_weights <- rep(rest_weight, rest_size)
  } else{
    rest_t <- NULL
    rest_weight <- 0
    rest_weights <- NULL
  }

  # interweave the samples so that they are in pairs in the resulting sequence of samples
  interweave <- c(seq(from = 1, to = sample_size-1, length.out = sample_size/2),
                  seq(from = 2, to = sample_size, length.out = sample_size/2))
  sample_t <- rep(NA, sample_size)
  sample_t[interweave] <- c(rej_t, rest_t)
  sample_weights <- rep(NA, sample_size)
  sample_weights[interweave] <- c(rej_weights, rest_weights)

  # convert the samples of T_j to response vector y
  sample_vjy <- tj_to_vjy(sample_t, y_Pi_Xnoj_res_norm2, df)
  y_results <- vjy_to_y(sample_vjy, j, y.pack, X.pack)

  return(list(y_samples = y_results$y_samples, sample_weights = sample_weights,
              weights = list(rej_weight = rej_weight, rest_weight = rest_weight),
              sigmahat_X_res = y_results$sigmahat_X_res,
              sigmahat_XXk_res = y_results$sigmahat_XXk_res))
}

# the set of tval[j] where Hj marginally rejected conditional on Sj
where_pval_rej <- function(j, y.pack){
  sigmahat <- sqrt(y.pack$y_Pi_X_res_norm2 / y.pack$df)
  tj_obs <- as.numeric(y.pack$vjy_obs[j] / sigmahat)

  tval <- abs(tj_obs)

  rej_set <- list(left = c(-Inf, tval), right = c(-tval, Inf))

  return(rej_set)
}

# approximately find the the set of tval[j] (and vj*y) where knockoff reject j conditional on Sj
where_kn_rej <- function(kn_alpha, j, y.pack, X.pack,
                         statistic,
                         method = "local_lin_reg"){
  df <- y.pack$df
  y_Pi_Xnoj_res_norm2 <- y.pack$y_Pi_Xnoj_res_norm2[j]

  # find the approximated knockoff rejection region in vjy
  kn_stat_samples <- kn_stat_sampling(kn_alpha, j, y.pack, X.pack, node_num = 10,
                                      statistic)

  vjy_region <- region_F1geqF2(kn_stat_samples$vjy_nodes,
                               kn_stat_samples$kn_abs_stat_j,
                               kn_stat_samples$kn_stat_thr,
                               method)
  left_vjy <- vjy_region$left
  right_vjy <- vjy_region$right

  # reorganize the results
  if(length(left_vjy)>0){
    # concatenate continuous intervals
    left <- left_vjy[1]
    right <- NULL
    if(length(left_vjy) >= 2){
      for(i in 1:(length(left_vjy)-1)){
        if(right_vjy[i] != left_vjy[i+1]){
          right <- c(right, right_vjy[i])
          left <- c(left, left_vjy[i+1])
        }
      }
    }
    right <- c(right, tail(right_vjy, 1))
    left_vjy <- left
    right_vjy <- right

    # convert vjy to t-val
    left_t <- vjy_to_tj(left_vjy, y_Pi_Xnoj_res_norm2, df)
    right_t <- vjy_to_tj(right_vjy, y_Pi_Xnoj_res_norm2, df)

  } else{
    left_vjy <- NULL
    right_vjy <- NULL
    left_t <- NULL
    right_t <- NULL
  }

  rej_set <- list(left = left_t, right = right_t,
                  left_vjy = left_vjy, right_vjy = right_vjy)
}

# generating nodes on vj*y. Compute the kn-stat Wj and the kn-rejection lower
# bound W_{(\hat k)} on the them.
kn_stat_sampling <- function(kn_alpha, j, y.pack, X.pack, node_num = 10,
                             statistic){
  df <- y.pack$df
  y_Pi_Xnoj_res_norm2 <- y.pack$y_Pi_Xnoj_res_norm2[j]
  vjy_obs <- y.pack$vjy_obs[j]

  # compute the bounds of the finite interval to interpolate/LLR
  tj_bound <- qt(1e-3, df = df, lower.tail = F)
  vjy_bound <- max(abs(tj_to_vjy(tj_bound, y_Pi_Xnoj_res_norm2, df)), abs(vjy_obs))
  vjy_nodes <- seq(-vjy_bound, vjy_bound, length.out = node_num)

  y_results <- vjy_to_y(vjy_nodes, j, y.pack, X.pack)
  y_nodes <- y_results$y_samples

  # compute the kn-stat Wj, the kn-rejection lower bound W_{(hat k)} on the nodes
  kn_abs_stat_j <- rep(NA, node_num)
  kn_stat_thrs <- rep(NA, node_num)
  for(node_i in 1:node_num){
    if("sigma_tilde" %in% names(formals(statistic))){
      kn_stats <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, node_i],
                            sigma_tilde = y_results$sigmahat_XXk_res[node_i])
    } else{
      kn_stats <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, node_i])
    }

    kn_abs_stat_j[node_i] <- abs(kn_stats[j])
    kn_stat_thrs[node_i] <- (kn.select(kn_stats, kn_alpha,
                                       selective = T, early_stop = 1))$W_k_hat
  }

  # return the values
  return(list(vjy_nodes = vjy_nodes, kn_abs_stat_j = kn_abs_stat_j,
              kn_stat_thrs = kn_stat_thrs))
}


region_F1geqF2 <- function(x, y1, y2, method){
  # linear interpolation
  if(method == "lin_interp"){
    int_num = length(x) - 1 # number of intervals given by the interpolation
    left <- NULL # left bounds of intervals of vj*y where F1 >= F2
    right <- NULL # right bounds ...

    # in each interval between consecutive x nodes, find the region where
    # linearly interpolate F1 >= minimal of F2 at the end points
    for(int_i in 1:int_num){
      # set F2 to be the minimum at the two boundary nodes
      thr <- min(y2[int_i], y2[int_i+1])

      # solve the linear inequality y2 >= y1
      if(y1[int_i] >= thr & y1[int_i+1] >= thr){
        left <- c(left, x[int_i])
        right <- c(right, x[int_i+1])
      } else if(y1[int_i] > thr | y1[int_i+1] > thr){
        intersect <- (thr - y1[int_i]) / (y1[int_i+1] - y1[int_i]) *
          (x[int_i+1] - x[int_i]) + x[int_i]
        if(y1[int_i] < y1[int_i+1]){
          left <- c(left, intersect)
          right <- c(right, x[int_i+1])
        } else{
          left <- c(left, x[int_i])
          right <- c(right, intersect)
        }
      }
    }
  }
  # local linear regression
  else if(method == "local_lin_reg"){
    bandwidth <- max(abs(diff(x)))    # set the band width as the distance between nodes
    x_range <- c(min(x), max(x))
    n_eval <- 100

    # local linear regression with Gaussian kernel
    F1 <- KernSmooth::locpoly(x, y1, bandwidth = bandwidth, degree = 1, gridsize = n_eval, range.x = x_range)
    F2 <- KernSmooth::locpoly(x, y2, bandwidth = bandwidth*2, degree = 1, gridsize = n_eval, range.x = x_range)

    # find the desired region
    F1geqF2 <- F1$y >= F2$y
    # extend the region by one evaluation node
    F1geqF2 <- (F1geqF2 | c(F, F1geqF2[1:(n_eval-1)]) | c(F1geqF2[2:n_eval], F))
    # identify the boundary
    x_eval <- F1$x
    left <- x_eval[which(F1geqF2 & !c(F, F1geqF2[1:(n_eval-1)]))]
    right <- x_eval[which(F1geqF2 & !c(F1geqF2[2:n_eval], F))]
  }

  return(list(left = left, right = right))
}


vjy_to_tj <- function(vjy, y_Pi_Xnoj_res_norm2, df){
  tj <- vjy / sqrt((y_Pi_Xnoj_res_norm2 - vjy^2) / df)
}

tj_to_vjy <- function(tj, y_Pi_Xnoj_res_norm2, df){
  vjy <-  tj * sqrt(y_Pi_Xnoj_res_norm2 / (tj^2 + df))
}

# recover y from vjy conditional on Sj
vjy_to_y <- function(sample_vjy, j, y.pack, X.pack){
  # retrieve data about X
  Xk_dim <- y.pack$p
  XXk_res_dim <- y.pack$n - 2*y.pack$p
  vj <- X.pack$vj

  # retrieve data about the observed y
  y_Pi_Xnoj_res_norm2 <- y.pack$y_Pi_Xnoj_res_norm2[j]
  y_Pi_Xnoj <- y.pack$y_Pi_Xnoj[, j]

  sample_size <- length(sample_vjy)

  # square norm of y in residue of X
  y_Pi_X_res_norm2_sample <- y_Pi_Xnoj_res_norm2 - sample_vjy^2

  ## sample the y projected on to the residue of X
  # coordinates of y in (Pi_{X}^\perp Xk) under the QR orthogonal basis
  y_Pi_X_res_Xk_sample <- matrix(rnorm(Xk_dim * sample_size), ncol = sample_size)
  # square norm of y in (Pi_{[X Xk]}^\perp), note the direction is ancillary
  y_Pi_XXk_res_norm2_sample <- rchisq(sample_size, df = XXk_res_dim)

  # scale the above two variables to satisfy |Pi_{X}^\perp y|
  y_res_X_norm2_scale <- (colSums(y_Pi_X_res_Xk_sample^2) + y_Pi_XXk_res_norm2_sample) / y_Pi_X_res_norm2_sample
  y_Pi_X_res_Xk_sample <- scale(y_Pi_X_res_Xk_sample, center = F, scale = sqrt(y_res_X_norm2_scale))
  y_Pi_XXk_res_norm2_sample <- y_Pi_XXk_res_norm2_sample / y_res_X_norm2_scale

  # assemble y in (Pi_{X}^\perp), direction in (Pi_{[X Xk]}^\perp) is picked arbitrarily as XXk_res_unit
  y_Pi_X_res_sample <- X.pack$X_res_Xk_basis %*% y_Pi_X_res_Xk_sample +
    X.pack$XXk_res_unit %*% matrix(sqrt(y_Pi_XXk_res_norm2_sample), nrow = 1)

  # assemble y
  y_samples <- matrix(rep(y_Pi_Xnoj, sample_size), ncol = sample_size) + vj %*% matrix(sample_vjy, nrow = 1) + y_Pi_X_res_sample

  # get the sigmahat for free, used for computing t-value and statistic
  sigmahat_X_res <- sqrt(y_Pi_X_res_norm2_sample / (y.pack$n - y.pack$p))
  sigmahat_XXk_res <- sqrt(y_Pi_XXk_res_norm2_sample / (y.pack$n - 2*y.pack$p))

  return(list(y_samples = y_samples,
              sigmahat_X_res = sigmahat_X_res, sigmahat_XXk_res = sigmahat_XXk_res))
}

#' @importFrom utils tail
# generate samples uniformly distribution on a union of intervals (lowers, uppers)
runif_intervals <- function(lowers, uppers, sample_size){
  itv_thres <- cumsum(uppers - lowers)
  total_mass <- tail(itv_thres, 1)

  samples <- runif(sample_size, min = 0, max = total_mass)
  samples <- sapply(samples, function(sample){
    ind <- min(which(sample <= itv_thres))
    sample <- sample - c(0, itv_thres)[ind] + lowers[ind]
  })

  return(samples)
}






