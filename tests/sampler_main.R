# call cKnockoff_j() multiple times and aggregate the results
invest_ineq <- function(n, p, # problem size
                        mc_size, # Monte-Carlo sample size
                        X_type, # type of design matrix
                        intercept, # whether fit with intercept
                        target, # power of BH at level 0.2, this determines the signal strength
                        alpha, # FDR level
                        alt_num, # number of non-nulls
                        y_seed, # random seed for reproducibility
                        invest_j, # feature j to invest
                        naive_sampler = T, # shoule we use the naive sampler or our actual sampler
                        sample_coupling = F, # should we couple the Monte-Carlo samples
                        kn_reg_method = "local_lin_reg", # how do we fit the |W_j| curve, alternative is "lin_interp", linear interpolation
                        statistic = stat.glmnet_coefdiff_lm, # feature statistics
                        noise = quote(rnorm(n)), # noise in the linear model
                        n_cores = 10 # number of cores in parallel computing
){
  X_seed <- 1
  X <- gene_X(X_type, n, p, X_seed)

  if(X_type == "X_block"){
    pi1 <- 10 / p
    mu1 <- BH_lm_calib(X, pi1, noise = quote(rnorm(n)), "fix", 1, side = "two", nreps = 200, alpha = alpha,
                       target = target, n_cores = n_cores)

    beta <- genmu(p, pi1, mu1, "fix", 1)
  } else{
    pi1 <- alt_num / p
    mu1 <- BH_lm_calib(X, pi1, noise = quote(rnorm(n)), "random", 1, side = "two", nreps = 200, alpha = alpha,
                       target = target, n_cores = n_cores)

    beta <- genmu(p, pi1, mu1, "random", 1)
  }
  H0 <- beta == 0

  # registerDoParallel(n_cores)

  # precompute the matrices related to X
  X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))
  X.pack <- process_X(X, intercept = intercept)


  set.seed(y_seed)
  y <- X %*% beta + eval(noise)
  # browser()

  ineq_invest <- cKnockoff_j(y, X, X.pack,
                             intercept = intercept,
                             alpha = alpha,
                             H0, seed = y_seed+100,
                             mc_size = mc_size,
                             invest_j = invest_j,
                             naive_sampler = naive_sampler,
                             kn_reg_method = kn_reg_method,
                             statistic = statistic,
                             sample_coupling = sample_coupling)


  rejected <- sum((ineq_invest$p_ineq$ineq_L - ineq_invest$p_ineq$ineq_R) *
                    ineq_invest$samples.weight) <= 0

  z_norm2 <- sum(y^2) - sum((lm(y ~ X[, -invest_j])$fitted.values)^2)

  return(c(ineq_invest, list(X = X, y = y, z_norm2 = z_norm2,
                             alt = which(!H0), rejected = rejected)))
}


# cKnockoff for a particular feature j
cKnockoff_j <- function(y, X,
                        X.pack,
                        intercept = intercept,
                        alpha = 0.2,
                        H0, # a T/F vector indicating if H_j is true null
                        seed = 1,
                        mc_size = 100, # Monte-Carlo sample size
                        invest_j, # feature j to invest
                        naive_sampler = T, # shoule we use the naive sampler or our actual sampler
                        kn_reg_method = "local_lin_reg", # how do we fit the |W_j| curve, alternative is "lin_interp", linear interpolation
                        statistic, # feature statistics
                        sample_coupling = F # should we couple the Monte-Carlo samples
                        ){
  invest_mode <- T

  j <- invest_j

  X.org <- X
  y.org <- y

  y.data <- transform_y(X.pack, y, intercept = intercept, randomize = F)
  y.pack <- process_y(X.pack, y.data)

  X <- X.pack$X
  y <- y.pack$y.data$y
  df_X <- y.pack$y.data$df_X

  n <- NROW(X)
  p <- NCOL(X)

  # Xk_dir <- X.pack$X_res_Xk_basis[, 2]

  Xkn_full <- scale(cbind(X, X.pack$X_kn))

  if(naive_sampler){
    tvals_obs <- lm_to_t(y, X, sigmahat = sqrt(y.pack$RSS_X / df_X))$tvals
  } else{
    tvals_obs <- y_to_t(y, X.pack$vj_mat, sqrt(y.pack$RSS_X / df_X))
  }
  pvals_obs <- pvals_t(tvals_obs, df = df_X, side = "two") # marginal test qvals
  pval_obs <- pvals_obs[j]

  ineq_LR <- matrix(rep(0, 2*mc_size), 2)

  sigma_est <- sqrt(y.pack$RSS_Xnoj[j] / (df_X+1)) # overestimate sigma
  fit_on_rest <- glmnet::glmnet(X[, -j], y,
                                lambda = 2 * sigma_est / n,
                                intercept=T, standardize=F, standardize.response=F,
                                family="gaussian")
  y.pack$Xy_bias[j] <- t(X[, j]) %*% (X[, -j] %*% fit_on_rest$beta + rep(fit_on_rest$a0, each = n))
  cali_stats_obs <- abs(c(matrix(y, nrow = 1) %*% X) - y.pack$Xy_bias)

  # for investing
  if(invest_mode){
    p_ineq <- data.frame(pval = numeric(), tval = numeric(),
                         ineq_L = numeric(), ineq_R = numeric(),
                         kn_rej = logical(), nokn_rej = logical(),
                         kept = logical())
    pval_mat <- matrix(NA, p, mc_size)
    kn_stats_mc <- matrix(NA, mc_size, p)
    eskn_rej_mc <- matrix(NA, mc_size, p)
    vjy_mc <- rep(NA, mc_size)
    # SRL_est_mc <- rep(NA, mc_size)
    kn_stat_thres <- rep(NA, mc_size)
    # pi0_mc <- rep(NA, mc_size)
    # Rkn_mc <- rep(NA, mc_size)
    # Reskn_mc <- rep(NA, mc_size)
    # qval_invest <- pval_obs
    y_Pi_Xk1_mc <- rep(NA, mc_size)
    pval_mc <- rep(NA, mc_size)
    weights_mc <- NULL
  }


  # knockoff and dBY rejection
  kn_stats_obs <- statistic(X, X.pack$X_kn, y,
                            sigma_tilde = sqrt(y.pack$y.data$RSS_XXk
                                               / y.pack$y.data$df_XXk))
  kn_selected <- kn.select(kn_stats_obs, alpha,
                           selective = T, early_stop = 0)$selected

  init_selected <- j %in% kn_selected


  ## Monte Carlo sample for RB
  # # naive sampler
  if(naive_sampler){
    samples.res <- y_condj_sample(y.org, X.org, j, mc_size, seed = seed+j)
    samples.y <- sapply(1:mc_size, function(mc_i){
      as.vector(t(samples.res[, mc_i]) %*% X.pack$XXk.org.basis)
    })
    samples.weight <- rep(1, mc_size)
    sigmahat.X_res <- sqrt((colSums(samples.res^2) - colSums(samples.y[1:p, ]^2)) /
                             y.pack$y.data$df_X)
    sigmahat.XXk_res <- sqrt((colSums(samples.res^2) - colSums(samples.y^2)) /
                               y.pack$y.data$df_XXk)
  }

  # rejection region of the new sampler
  kn_stat_nodes <- kn_stat_sampling(alpha, j, y.pack, X.pack,
                                    node_num = 10,
                                    statistic)
  kn_rej_set <- where_kn_rej(alpha, j, y.pack,
                             X.pack = NULL,
                             statistic, method = kn_reg_method,
                             kn_stat_samples = kn_stat_nodes)


  cali_rej_set <- where_cali_rej(j, y.pack, X.pack)
  sample_region <- interval_union(kn_rej_set, cali_rej_set)

  samples.res <- y_sampler_cond_Sj(mc_size, cali_rej_set, kn_rej_set, j,
                                  y.pack, X.pack, sample_coupling)
  if(!naive_sampler){
    samples.y <- samples.res$samples.y
    samples.weight <- samples.res$samples.weight
    sigmahat.X_res <- sqrt(samples.res$samples.RSS_X / samples.res$df_X)
    sigmahat.XXk_res <- sqrt(samples.res$samples.RSS_XXk / samples.res$df_XXk)

    ineq_bounds <- get_Ej_bound(alpha, p, samples.res$weights, sample_coupling)
  }

  if(invest_mode){
    weights_mc <- samples.weight
  }

  Ej_mc <- rep(NA, mc_size)

  for(mc_i in 1:mc_size){

    # compute knockoff stat
    kn_stats <- statistic(X, X.pack$X_kn, samples.y[, mc_i], sigmahat.XXk_res[mc_i])

    ## compute weights
    # make rejection and fdp estimation based on knockoff statistics
    # selective SeqStep: selective = TRUE; SeqStep: selective = FALSE
    kn_result <- kn.select(kn_stats, alpha,
                           selective = T, early_stop = 1)
    kn_selected_mc <- kn_result$selected

    # compute b_j
    kn_null_num <- (kn_result$fdp_est * length(kn_selected_mc)) %>% max(1)
    b_j <- alpha * (j %in% kn_selected_mc) / kn_null_num

    # for investing
    if(invest_mode){
      eskn_rej_mc[mc_i, ] <- (1:p) %in% kn_selected_mc
      kn_stat_thres[mc_i] <- kn_result$W_k_hat
      # pi0_mc[mc_i] <- kn_result$fdp_est
      # Reskn_mc[mc_i] <- length(kn_selected_mc)
    }

    ## compute DP_j
    # p-values
    if(naive_sampler){
      tvals_mc <- lm_to_t(samples.y[, mc_i], X, sigmahat = sigmahat.X_res[mc_i])$tvals
    } else{
      tvals_mc <- y_to_t(samples.y[, mc_i], X.pack$vj_mat, sigmahat.X_res[mc_i])
    }
    pvals_mc <- pvals_t(tvals_mc, df_X, side = "two")

    # test if t-value is computed correctly
    if(!naive_sampler && F){
      t_naive <- lm_to_t(samples.y[, mc_i], X, sigmahat = sigmahat.X_res[mc_i])$tvals
      t_new <- y_to_t(samples.y[, mc_i], X.pack$vj_mat, sigmahat.X_res[mc_i])
      print(max(abs(t_new - t_naive)) < 1e-12)

      sigma_tilde_naive <- sqrt((sum(samples.y[, mc_i]^2) -
                                   sum((lm(samples.y[, mc_i] ~ cbind(X.pack$X, X.pack$X_kn))$fitted.values)^2)) / (n - 2*p))
      sigma_tilde_new <- sigmahat.XXk_res[mc_i]
      print(max(sigma_tilde_new - sigma_tilde_naive))
      browser()
    }

    # knockoff rejections
    kn_result <- kn.select(kn_stats, alpha,
                           selective = T, early_stop = 0)
    kn_selected_mc <- kn_result$selected

    # # for investing
    # if(invest_mode){
    #     Rkn_mc[mc_i] <- length(kn_selected_mc)
    # }

    # R hat
    selected_mc <- kn_selected_mc

    # compute calibration selection
    cali_stat_mc <- abs(sum(X[, j] * samples.y[, mc_i]) - y.pack$Xy_bias[j])
    cali_selected_mc <- ifelse(cali_stat_mc >= cali_stats_obs[j], j, 0)

    # compute DP_j
    DP_j <- (j %in% union(kn_selected_mc, cali_selected_mc)) / length(union(selected_mc, j))

    ineq_LR[, mc_i] <- c(DP_j, b_j) * samples.weight[mc_i]
    # for investing
    if(invest_mode){

      pval_invest <- pvals_t(tvals_mc, df_X, side = "right")
      # if(pval_invest[j] > 0.99) browser()

      p_ineq <- add_row(p_ineq,
                        pval = pval_invest[j],
                        tval = tvals_mc[j],
                        ineq_L = DP_j,
                        ineq_R = b_j,
                        kn_rej = (j %in% kn_selected_mc),
                        nokn_rej = (j %in% cali_selected_mc),
                        kept = x_in_intervals(tvals_mc[j], sample_region))

      pval_mat[, mc_i] <- pval_invest
      kn_stats_mc[mc_i, ] <- kn_stats
      vjy_mc[mc_i] <- t(X.pack$vj_mat[, j]) %*% samples.y[, mc_i]
      # SRL_est_mc[mc_i] <- (max(abs(t(Xkn_full) %*% samples.y[, mc_i])) +
      #                          max(abs(t(Xkn_full[, j] %*% samples.y[, mc_i])),
      #                              abs(t(Xkn_full[, j + p] %*% samples.y[, mc_i])))) / 2 / sqrt(n-1)

      y_Pi_Xk1_mc[mc_i] <- samples.y[p+1, mc_i]
      pval_mc[mc_i] <- pval_invest[j]
    }

  }

  # test the correctness of the importance sampler and weights by comparing with the naive sampler
  inds_in_rej <- sapply(p_ineq$tval, function(x){
    x_in_intervals(x, sample_region)
  })
  inds_in_rej <- which(inds_in_rej)
  num_in_rej <- length(inds_in_rej)
  ineq_in_rej <- (p_ineq$ineq_L[inds_in_rej] - p_ineq$ineq_R[inds_in_rej]) * samples.weight[inds_in_rej] * num_in_rej / mc_size
  print(paste0("sample_mass (empirical): ", num_in_rej / mc_size))
  print(paste0("Ej: ", mean(ineq_in_rej)))
  if(!naive_sampler){
    if(sample_coupling){
      ineq_combine <- (ineq_in_rej[seq(1, mc_size-1, length.out = mc_size/2)] + ineq_in_rej[seq(2, mc_size, length.out = mc_size/2)]) / 2
    } else{
      ineq_combine <- ineq_in_rej
    }
    ineq_bounds <- get_Ej_bound(alpha, p, samples.res$weights, sample_coupling)
    print(paste0("outside seq bound: ", sum((ineq_combine > ineq_bounds$upper) | (ineq_combine < ineq_bounds$lower))))
  }

  ineq_LR_mean <- rowMeans(ineq_LR)
  ineq_selected <- (ineq_LR_mean[1] <= ineq_LR_mean[2])


  return(list(p_ineq = p_ineq, pval_mat = pval_mat, kn_stats_mc = kn_stats_mc,
              eskn_rej_mc = eskn_rej_mc, vjy_mc = vjy_mc,
              # SRL_est_mc = SRL_est_mc,
              kn_stat_thres = kn_stat_thres, samples.weight = samples.weight,
              # qval_invest = qval_invest,
              # pi0_mc = pi0_mc, Rkn_mc = Rkn_mc, Reskn_mc = Reskn_mc,
              # init_selected = init_selected, ineq_selected = ineq_selected,
              kn_rej_set = kn_rej_set, cali_rej_set = cali_rej_set,
              sample_region = sample_region,
              kn_stat_nodes = kn_stat_nodes,
              y_Pi_Xk1_mc = y_Pi_Xk1_mc, pval_mc = pval_mc,
              weights_mc = weights_mc))
}



