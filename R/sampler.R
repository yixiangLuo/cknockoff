
# sampling from the p-value (calibration) rejection set union the knockoff rejection set
y_sampler_cond_Sj <- function(sample_size, cali_rej_reg, kn_rej_reg, j,
                              y.pack, X.pack, sample_coupling){
  df_X <- y.pack$y.data$df_X
  RSS_Xnoj <- y.pack$RSS_Xnoj[j]

  # mass of calibration rejection set
  rej_set <- cali_rej_reg
  if(!is.null(rej_set$left)){
    cali_rej_lb <- pt(rej_set$left, df = df_X, lower.tail = T)
    cali_rej_ub <- pt(rej_set$right, df = df_X, lower.tail = T)
    rej_mass <- sum(cali_rej_ub - cali_rej_lb)
  } else{
    rej_mass <- 0
  }

  # mass of (kn rejection set) - (calibration rejection set)
  rest_set <- interval_minus(kn_rej_reg, cali_rej_reg)
  if(!is.null(rest_set$left)){
    cali_rest_lb <- pt(rest_set$left, df = df_X, lower.tail = T)
    cali_rest_ub <- pt(rest_set$right, df = df_X, lower.tail = T)
    rest_mass <- sum(cali_rest_ub - cali_rest_lb)
  } else{
    rest_mass <- 0
  }

  total_mass <- rej_mass + rest_mass

  # can decide rejecting H_j or not based on the mass of the two sets
  if(total_mass == 0){
    return(0)
  }
  if(sample_coupling && rest_mass == 0){
    # sample_coupling <- F
    return(1)
  }
  if(rej_mass == 0){
    return(-1)
  }

  if(sample_coupling){
    rej_sample_prop <- 1/2
  } else{
    rej_sample_prop <- min(1, ceiling((rej_mass/total_mass * sample_size)) / sample_size)
  }

  # sampling in the calibration rejection set and compute the
  # importance sampling weights
  rej_size <- round(sample_size * rej_sample_prop)
  if(rej_size > 0){
    rej_p <- runif_intervals(cali_rej_lb, cali_rej_ub, rej_size)
    rej_t <- qt(rej_p, df = df_X, lower.tail = T)

    rej_weight <- rej_mass / rej_sample_prop
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
    rest_p <- runif_intervals(cali_rest_lb, cali_rest_ub, rest_size)
    rest_t <- qt(rest_p, df = df_X, lower.tail = T)

    rest_weight <- rest_mass / (1 - rej_sample_prop)
    rest_weights <- rep(rest_weight, rest_size)
  } else{
    rest_t <- NULL
    rest_weight <- 0
    rest_weights <- NULL
  }

  if(!sample_coupling){
    # randomly shuffle the samples
    if(sample_size > 1){
      shuffle <- c(1, sample(2:sample_size)) # make sure at least one sample from rej_set is seen
      sample_t <- c(rej_t, rest_t)[shuffle]
      sample_weights <- c(rej_weights, rest_weights)[shuffle]
    }
  } else{
    # interweave the samples so that they are in pairs in the resulting sequence of samples
    interweave <- c(seq(from = 1, to = sample_size-1, length.out = sample_size/2),
                    seq(from = 2, to = sample_size, length.out = sample_size/2))
    sample_t <- rep(NA, sample_size)
    sample_t[interweave] <- c(rej_t, rest_t)
    sample_weights <- rep(NA, sample_size)
    sample_weights[interweave] <- c(rej_weights, rest_weights)
  }

  # convert the samples of T_j to response vector y
  sample_vjy <- tj_to_vjy(sample_t, RSS_Xnoj, df_X)
  y_results <- vjy_to_y(sample_vjy, j, y.pack, X.pack)

  return(list(samples.y = y_results$y,
              samples.RSS_X = y_results$RSS_X,
              samples.RSS_XXk = y_results$RSS_XXk,
              df_X = df_X,
              df_XXk = y_results$df_XXk,
              samples.weight = sample_weights,
              weights = list(rej_weight = rej_weight, rest_weight = rest_weight)))
}


# determine should we couple the Monte-Carlo samples in computing E_j to improve
# efficiency
couple_samples <- function(df_X, cali_rej_set, kn_rej_set){

  # compute mass of the calibration rejection set
  rej_set <- cali_rej_set
  if(!is.null(rej_set$left)){
    cali_rej_lb <- pt(rej_set$left, df = df_X, lower.tail = T)
    cali_rej_ub <- pt(rej_set$right, df = df_X, lower.tail = T)
    rej_mass <- sum(cali_rej_ub - cali_rej_lb)
  } else{
    return(F)
  }

  # compute mass of the rest rejection set
  rest_set <- interval_minus(kn_rej_set, cali_rej_set)
  if(!is.null(rest_set$left)){
    cali_rest_lb <- pt(rest_set$left, df = df_X, lower.tail = T)
    cali_rest_ub <- pt(rest_set$right, df = df_X, lower.tail = T)
    rest_mass <- sum(cali_rest_ub - cali_rest_lb)
  } else{
    return(F)
  }

  # if rej_mass and rest_mass are not far from each other, couple the samples to
  # reduce vairance, otherwise don't to avoid over sampling in a relatively
  # small region
  threshold <- 3
  sample_coupling <- (rej_mass / rest_mass <= threshold) &&
    (rest_mass / rej_mass <= threshold)

  return(sample_coupling)
}

# find the set of tval[j] where Hj rejected by calibration conditional on Sj
where_cali_rej <- function(j, y.pack, X.pack){
  df_X <- y.pack$y.data$df_X
  RSS_Xnoj <- y.pack$RSS_Xnoj[j]

  Xjy_obs <- sum(X.pack$X[, j] * y.pack$y.data$y)
  Xjy_obs_counter <- 2 * y.pack$Xy_bias[j] - Xjy_obs
  X_proj_v <- sum(X.pack$X[, j] * X.pack$vj_mat[, j])
  vjy_obs_counter <- (Xjy_obs_counter - (Xjy_obs - y.pack$vjy_obs[j] * X_proj_v)) / X_proj_v

  vjy_rej_points <- sort(c(y.pack$vjy_obs[j], vjy_obs_counter))

  # max boundary of the values
  t_bound <- abs(qt(1e-14, df = df_X))
  vjy_bound <- tj_to_vjy(t_bound, RSS_Xnoj, df_X)

  rej_set <- list(left = NULL, right = NULL)

  # convert vjy to t-val
  if(vjy_rej_points[1] > -vjy_bound){
    rej_set$left <- c(rej_set$left, -t_bound)
    rej_set$right <- c(rej_set$right, vjy_to_tj(vjy_rej_points[1], RSS_Xnoj, df_X))
  }
  if(vjy_rej_points[2] < vjy_bound){
    rej_set$left <- c(rej_set$left, vjy_to_tj(vjy_rej_points[2], RSS_Xnoj, df_X))
    rej_set$right <- c(rej_set$right, t_bound)
  }

  return(rej_set)
}

# approximately find the the set of tval[j] (and vj*y) where knockoff reject j conditional on Sj
where_kn_rej <- function(kn_alpha, j, y.pack, X.pack,
                         statistic,
                         method = "local_lin_reg",
                         kn_stat_samples = NULL){
  df_X <- y.pack$y.data$df_X
  RSS_Xnoj <- y.pack$RSS_Xnoj[j]

  # find the approximated knockoff rejection region in vjy
  if(is.null(kn_stat_samples)){
    kn_stat_samples <- kn_stat_sampling(kn_alpha, j, y.pack, X.pack, node_num = 10,
                                        statistic)
  }

  # solve for the region where kn_abs_stat_j >= kn_stat_thr
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
    left_t <- vjy_to_tj(left_vjy, RSS_Xnoj, df_X)
    right_t <- vjy_to_tj(right_vjy, RSS_Xnoj, df_X)


    t_bound <- abs(qt(1e-14, df = df_X))
    boundary_t <- vjy_to_tj(c(kn_stat_samples$vjy_nodes[1],
                              kn_stat_samples$vjy_nodes[length(kn_stat_samples$vjy_nodes)]),
                            RSS_Xnoj, df_X)
    if(left_t[1] <= min(boundary_t)+1e-7){
      left_t[1] <- min(left_t[1], -t_bound)
    }
    if(right_t[length(right_t)] >= max(boundary_t)-1e-7){
      right_t[length(right_t)] <- max(right_t[length(right_t)], t_bound)
    }


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
  df_X <- y.pack$y.data$df_X
  RSS_Xnoj <- y.pack$RSS_Xnoj[j]
  vjy_obs <- y.pack$vjy_obs[j]

  # compute the bounds of the finite interval to interpolate/LLR
  tj_bound <- qt(5e-3, df = df_X, lower.tail = F)
  vjy_bound <- max(abs(tj_to_vjy(tj_bound, RSS_Xnoj, df_X)), abs(vjy_obs))
  vjy_nodes <- seq(-vjy_bound, vjy_bound, length.out = node_num)

  y_results <- vjy_to_y(vjy_nodes, j, y.pack, X.pack)
  y_nodes <- y_results$y
  sigmahat.XXk_res <- sqrt(y_results$RSS_XXk / y_results$df_XXk)

  # compute the kn-stat Wj, the kn-rejection lower bound W_{(hat k)} on the nodes
  kn_abs_stat_j <- rep(NA, node_num)
  kn_stat_thrs <- rep(NA, node_num)
  for(node_i in 1:node_num){
    if("sigma_tilde" %in% names(formals(statistic))){
      kn_stats <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, node_i],
                            sigma_tilde = sigmahat.XXk_res[node_i])
    } else{
      kn_stats <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, node_i])
    }

    kn_abs_stat_j[node_i] <- abs(kn_stats[j])
    kn_stat_thrs[node_i] <- (kn.select(kn_stats, kn_alpha,
                                       selective = T, early_stop = T))$W_k_hat
  }

  # return the values
  return(list(vjy_nodes = vjy_nodes, kn_abs_stat_j = kn_abs_stat_j,
              kn_stat_thrs = kn_stat_thrs))
}

# solve for the region where F1 >= F2
region_F1geqF2 <- function(x, y1, y2, method, n_eval = 100){
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
    n_eval <- n_eval

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


# convert vj^T y to corresponding t-statistic
vjy_to_tj <- function(vjy, RSS_Xnoj, df){
  tj <- vjy / sqrt((RSS_Xnoj - vjy^2) / df)
}

# convert t-statistic t_j to corresponding vj^T y
tj_to_vjy <- function(tj, RSS_Xnoj, df){
  vjy <-  tj * sqrt(RSS_Xnoj / (tj^2 + df))
}

# recover y from vjy conditional on Sj
vjy_to_y <- function(sample_vjy, j, y.pack, X.pack){
  # retrieve data about X
  p <- NCOL(X.pack$X_kn)
  df_XXk <- y.pack$y.data$df_XXk
  vj <- X.pack$vj_mat[, j]

  # retrieve data about the observed y
  RSS_Xnoj <- y.pack$RSS_Xnoj[j]
  y_Pi_Xnoj <- y.pack$y_Pi_Xnoj[, j]

  sample_size <- length(sample_vjy)

  # square norm of y in residue of X
  RSS_X_sample <- RSS_Xnoj - sample_vjy^2

  ## sample the y projected on to the residue of X
  # coordinates of y in (Pi_{X}^\perp Xk) in the QR orthonormal basis
  y_Pi_X_res_Xk_sample <- matrix(rnorm(p * sample_size), ncol = sample_size)
  # square norm of y in (Pi_{[X Xk]}^\perp), note the direction is ancillary
  RSS_XXk_sample <- rchisq(sample_size, df = df_XXk)

  # scale the above two variables to satisfy |Pi_{X}^\perp y|
  RSS_X_scale <- (colSums(y_Pi_X_res_Xk_sample^2) + RSS_XXk_sample) / RSS_X_sample
  y_Pi_X_res_Xk_sample <- scale(y_Pi_X_res_Xk_sample, center = F,
                                scale = sqrt(RSS_X_scale))
  RSS_XXk_sample <- RSS_XXk_sample / RSS_X_scale

  # assemble y in (Pi_{X}^\perp), direction in (Pi_{[X Xk]}^\perp) is suppressed
  y_Pi_X_res_sample <- rbind(matrix(0, nrow = p, ncol = sample_size),
                             y_Pi_X_res_Xk_sample)

  # assemble y
  y_samples <- matrix(rep(y_Pi_Xnoj, sample_size), ncol = sample_size) +
    vj %*% matrix(sample_vjy, nrow = 1) + y_Pi_X_res_sample

  return(list(y = y_samples,
              RSS_X = RSS_X_sample, RSS_XXk = RSS_XXk_sample,
              df_X = y.pack$y.data$df_X, df_XXk = df_XXk))
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






