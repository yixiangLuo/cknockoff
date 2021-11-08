#' @import stats methods
NULL


#' The cKnockoff procedure
#'
#' This function apply the cKnockoff procedure to the data under linear model,
#' selecting important features with FDR control. It has power dominating knockoff
#' via conditional calibration.
#'
#' @param X n-by-p matrix or data frame of predictors.
#' @param y response vector of length n.
#' @param knockoffs knockoffs method used to construct the knockoff matrix for X.
#' It should be a function taking a n-by-p matrix as input and returning a n-by-p
#' matrix of knockoff variables.
#' It is the same as the \code{knockoffs} parameter in \href{knockoff}{https://cran.r-project.org/web/packages/knockoff/}\code{::knockoff.filter}.
#' By default, \code{knockoff::create.fixed} is used.
#' @param statistic the knockoff W-statistics used to assess variable importance.
#' Any statistic function in \href{knockoff}{https://cran.r-project.org/web/packages/knockoff/}
#' can be passed via this parameter.
#' But please be aware of efficiency issue as this function will be called
#' repeatedly in cknockoff. We suggest use the statistic functions supplied by
#' our package.
#' By default, a lasso statistic from our package is used.
#' @param alpha target false discovery rate (default: 0.05).
#' @param kappa a real number between 0 and 1 controlling how should we hybridize
#' knockoff with Bonferroni in cKnockoff.
#' By default, kappa = 1, meaning we don't hybridize and use knockoff only.
#' @param n_cores the number of cores used in computing cKnockoff in parallel.
#' package \code{doParallel} is required if \code{n_cores} > 1.
#' Otherwise it's computed sequentially.
#' @param knockoff.type implies whether fixed-X or model-X knockoff is used.
#' Currently only supports fixed-X knockoff.
#' As a result, the parameter \code{knockoffs} has to be those for fixed-X.
#' @param prelim_result either a \code{knockoff.result} object returned by
#' \code{knockoff::knockoff.filter} or a \code{cknockoff.result} object returned
#' by \code{cknockoff}.
#' cknockoff can read the information from the knockoff result or previous
#' cknockoff result and make possibly more rejections under the same FDR control.
#' When supplied, all other parameters are not needed.
#' @param X.pack An object of class "cknockoff.X.pack" returned by \code{process_X}.
#' This is used only for simulations studies, where cknockoff is applied many
#' times to the same fixed X, to accelerate computation. General users should
#' ignore it and leave it as default = NULL.
#'
#' @return An object of class "cknockoff.result". This object is a list similar
#' to the "knockoff.result" object, containing essentially the same information:
#'  \item{X}{matrix of original variables (scaled and possibly augmented)}
#'  \item{Xk}{matrix of knockoff variables (cooresponding the the returned X)}
#'  \item{y}{response vector (possibly augmented)}
#'  \item{kn.statistic}{the knockoff W statistics}
#'  \item{selected}{named vector of selected variables}
#'  \item{record}{a list recording information used to assist computing from
#'   "prelim_result". Users should not worry or care about it.}
#'
#' @details
#'
#' #' The parameter \code{kappa} controls how Knockoff and Bonferroni is hybridized
#' in cKnockoff.
#' As kappa get closer to 1, cKnockoff behave more like Knockoff:
#' relying on the knockoff W statistics heavier.
#' While as kappa get closer to 0, cKnockoff behave more like Bonferroni:
#' relying on the p-values heavier.
#' In particular, cKnockoff(alpha; kappa) almost surely dominates
#' Knockoff(fdr = kappa \* alpha) union Bonferroni(fdr = (1-kappa) \* alpha).
#'
#' If parameter "prelim_result" is supplied, all the other parameters except
#' "n_cores" will be overwritten by the information extracted from it.
#'
#' If a \code{knockoff.result} object is supplied, cknockoff will return possibly
#' more selections on the same problem with the same FDR control.
#' To use, the function call of \code{knockoff::knockoff.filter} that generates
#' such \code{knockoff.result} object should specify parameters \code{statistic}
#' and \code{fdr} explicitly (rather than relying on the default).
#'
#' It's possible to even make more rejections based on a \code{cknockoff.result}
#' object. This is because cknockoff will only explore those promising features
#' and discards the others, say those with large p-values.
#' Calling cknockoff() on a \code{cknockoff.result} object would further
#' explore some most promising ones among the discarded features and (rarely)
#' possibly make several more rejections. You can do this recursively and FDR is
#' proved to be controlled whenever you decide to stop.
#' The computational time of each call should be similar.
#'
#' @references
#'
#' @export
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 3.5 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#' print(which(1:p %in% nonzero))
#'
#' # Basic usage
#' result <- cknockoff(X, y, alpha = 0.05, n_cores = 1)
#' print(result$selected)
#'
#' # improve knockoff and previous cknockoff result
#' library("knockoff")
#' kn.result <- knockoff.filter(X, y,
#'                              knockoffs = create.fixed,
#'                              statistic = stat.glmnet_coefdiff_lm,
#'                              fdr = 0.05)
#' print(kn.result$selected)
#' result <- cknockoff(prelim_result = kn.result, n_cores = 2)
#' print(result$selected)
#' result <- cknockoff(prelim_result = result)
#' print(result$selected)
cknockoff <- function(X, y,
                      knockoffs = knockoff::create.fixed,
                      statistic = stat.glmnet_coefdiff_lm,
                      alpha = 0.05,
                      kappa = 1,
                      n_cores = 1,
                      knockoff.type = c("fixed", "model"),
                      prelim_result = NULL,
                      X.pack = NULL,
                      Rhat_level = 0,
                      cali_with_kn = F,
                      cali_var = "pval"){

  mc_rounds <- 5
  mc_size <- 100

  Rhat_max_num <- 5
  Rhat_calc_max_step <- 4

  # knockoff.type <- match.arg(knockoff.type)

  # cknockoff.call <- match.call.defaults()
  # cknockoff.call[[1]] <- as.symbol("parse_args")
  # args <- eval(cknockoff.call)
  envir <- parent.frame()
  args <- parse_args(X, y, knockoffs, statistic,
                     alpha, kappa, n_cores, knockoff.type,
                     prelim_result, X.pack,
                     envir = envir)

  check_args(args)

  args <- process_args(args)
  X.pack <- args$X.pack
  y.pack <- args$y.pack
  X <- X.pack$X
  y <- y.pack$y
  statistic <- args$statistic
  alpha <- args$alpha
  kappa <- args$kappa
  n_cores <- args$n_cores
  record <- args$record

  forall <- args$parallel$iterator
  `%exec%` <- args$parallel$connector

  rm(args)

  if(!requireNamespace("KernSmooth", quietly=T)) {
    warning("KernSmooth is not installed. \n",
            "We use local linear regression from KernSmooth to assist Monte-Carlo sampling. \n",
            "Linear interpolation will be used instead. \n",
            "It might make the results a bit more sensitive to random noise.", call. = F, immediate. = T)
    kn_region_method <- "lin_interp"
  } else{
    kn_region_method <- "local_lin_reg"
  }

  n <- NROW(X)
  p <- NCOL(X)
  df <- n - p

  tvals_obs <- y_to_t(y, X.pack$vj_mat, sqrt(y.pack$y_Pi_X_res_norm2 / df))
  pvals_obs <- pvals_t(tvals_obs, df = df, side = "two")


  # knockoff and Bonferroni rejection
  if(is.null(record)){
    if("sigma_tilde" %in% names(formals(statistic))){
      kn_stats_obs <- statistic(X, X.pack$X_kn, y,
                                sigma_tilde = y.pack$sigmahat_XXk_res)
    } else{
      kn_stats_obs <- statistic(X, X.pack$X_kn, y)
    }

    kn_selected <- kn.select(kn_stats_obs, kappa * alpha,
                             selective = T, early_stop = F)$selected
  } else{
    kn_stats_obs <- record$kn.statistic
    kn_selected <- record$kn.selected
  }

  Bonf_selected <- which(pvals_obs <= (1-kappa) * alpha / p)
  init_selected <- union(kn_selected, Bonf_selected)

  if(!is.null(record) && record$iteration > 0){
    selected <- record$selected_so_far
  } else{
    selected <- init_selected
  }

  # find the promising features not selected or checked yet
  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha, record, selected, cali_with_kn = cali_with_kn)

  cali_stats_obs <- rep(NA, p)
  y_fit_noj_all <- matrix(NA, nrow = n, ncol = p)
  if(cali_var != "pval"){
    if(cali_var == "LM_lasso"){
      for(j in candidates){
        fit <- glmnet::glmnet(X[, -j], y, lambda=2*y.pack$sigmahat_XXk_res/n, intercept=T, standardize=F, standardize.response=F, family="gaussian")
        y_fit_noj_all[, j] <- as.vector(X[, -j] %*% fit$beta + rep(fit$a0, each = n))
      }
    }
    cali_stats_obs <- cali_statasitic(X.pack$X, X.pack$X_kn, y,
                                      statistic, sigma_tilde = y.pack$sigmahat_XXk_res,
                                      y_fit_noj_all, type = cali_var)
  }

  # check each hypothesis in sequence/parallel
  cali_selected <- forall(j = candidates, .options.multicore = list(preschedule = F)) %exec% {

    # prepare repeatedly used quantities
    pval_obs <- pvals_obs[j]
    mc_used <- 0
    ineq <- NULL
    ineq_combine <- NULL
    decision <- NA
    select_j <- F

    # prepare the sets where knockoff or a individual p-value test would reject j
    # conditional on Sj. This will be used in generating Monte-Carlo samples of y
    # conditional on Sj.
    kn_rej_set <- where_kn_rej(kappa * alpha, j, y.pack,
                               X.pack = list(X = X.pack$X, X_kn = X.pack$X_kn,
                                             vj = X.pack$vj_mat[, j],
                                             X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                             XXk_res_unit = X.pack$XXk_res_unit),
                               statistic, method = kn_region_method)
    if(cali_var == "pval"){
      pval_rej_set <- where_pval_rej(j, y.pack)
    } else{
      cali_thr <- cali_stats_obs[j]
      cali_rej_set <- where_cali_rej(j, y.pack,
                                     X.pack = list(X = X.pack$X, X_kn = X.pack$X_kn,
                                                   vj = X.pack$vj_mat[, j],
                                                   X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                                   XXk_res_unit = X.pack$XXk_res_unit),
                                     statistic,
                                     method = kn_region_method, cali_var,
                                     cali_thr = cali_thr,
                                     y_fit_noj_all = y_fit_noj_all)
    }


    if(Rhat_level > 0){
      relax_factor <- 2
      Rhat_to_rej_j_est <- Rhat_to_rej_j(alpha, kappa, n, p,
                                         pval_rej_set, kn_rej_set, cali_with_kn)
      Rhat_level_j <- ifelse(Rhat_to_rej_j_est <= Rhat_max_num * relax_factor &&
                               Rhat_to_rej_j_est >= 1 / relax_factor, 1, 0)
    } else{
      Rhat_level_j <- 0
    }
    Rhat_vec <- NULL


    # we divide the whole MC samples set into "mc_rounds" batches of size "mc_size"
    # and work on each batch sequentially for efficiency reason.
    for(mc_round in 1:mc_rounds){
      ineq <- c(ineq, rep(NA, mc_size))
      ineq_combine <- c(ineq_combine, rep(NA, mc_size/2))

      # generate Monte-Carlo samples of y conditional Sj for this batch
      if(cali_var == "pval"){
        sample_res <- y_sampler_cond_Sj(mc_size, pval_rej_set, kn_rej_set, j,
                                        y.pack,
                                        X.pack = list(vj = X.pack$vj_mat[, j],
                                                      X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                                      XXk_res_unit = X.pack$XXk_res_unit))
      } else{
        sample_res <- y_sampler_cond_Sj_newCali(mc_size, cali_rej_set, kn_rej_set, j,
                                                y.pack,
                                                X.pack = list(vj = X.pack$vj_mat[, j],
                                                              X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                                              XXk_res_unit = X.pack$XXk_res_unit))
      }


      # make decision if we can tell whether ineq <=0 based on the sampling region
      if(is.numeric(sample_res) && mc_round == 1){
        if(sample_res == 0) break
        else if(sample_res == 1){
          break
        }
        else if(sample_res == -1){
          select_j <- T
          break
        }
      }
      # retrive the samples
      y_cond <- sample_res$y_samples
      sample_weights <- sample_res$sample_weights

      sigmahat_X_res <- sample_res$sigmahat_X_res
      sigmahat_XXk_res <- sample_res$sigmahat_XXk_res

      # the upper and lower bound of ineq for constructing confidence sequence.
      ineq_bounds <- get_ineq_bound(kappa, alpha, p, sample_res$weights)

      # for each y MC sample in conditional calibration
      for(mc_i in 1:mc_size){

        # compute knockoff stat
        if("sigma_tilde" %in% names(formals(statistic))){
          kn_stat_mc <- statistic(X, X.pack$X_kn, y_cond[, mc_i],
                                  sigma_tilde = sigmahat_XXk_res[mc_i])
        } else{
          kn_stat_mc <- statistic(X, X.pack$X_kn, y_cond[, mc_i])
        }
        if(cali_var != "pval"){
          cali_stat_mc <- cali_statasitic(X, X.pack$X_kn, y_cond[, mc_i],
                                          statistic, sigma_tilde = sigmahat_XXk_res[mc_i],
                                          y_fit_noj_all, type = cali_var)[j]
        } else{
          cali_stat_mc <- NA
        }

        if(Rhat_level_j > 0 && length(Rhat_vec) >= 10){
          if(mean(1/Rhat_vec) >= min(1/1.2, 2/Rhat_to_rej_j_est) &&
             mean(ineq[1:mc_used]) > 0.08/mc_used){
            Rhat_level_j <- 0
          }
        }

        # if(Rhat_level_j > 0){
        #   X.pack_prepare <- X.pack
        # } else{
        #   X.pack_prepare <- NA
        # }
        fj_result <- calc_fj(j, y_cond[, mc_i], sample_weights[mc_i],
                             X.pack = X.pack, X.pack$vj_mat,
                             alpha, kappa, pval_obs,
                             statistic, kn_stat_mc, sigmahat_X_res[mc_i],
                             cali_with_kn, Rhat_level = Rhat_level_j,
                             Rhat_max_num, Rhat_calc_max_step,
                             cali_var, cali_stats_obs[j], cali_stat_mc)


        Rhat_vec <- c(Rhat_vec, fj_result$Rhat)

        # record the result
        mc_used <- mc_used + 1
        ineq[mc_used] <- fj_result$fj

        # compute the confidence interval and make decision
        decision <- make_decision(ineq[1:mc_used], ineq_bounds, threshold = 0,
                                  rej_alpha = min(0.05, alpha * max(1, length(init_selected)) / length(candidates)),
                                  accept_alpha = 0.05)
        if(decision$confident){ break }

      }

      # decide rejecting or not
      if(decision$confident){
        if(decision$reject){
          select_j <- T
        }
        break
      } else{
        if(!decision$reject){
          break
        } else{
          if(mc_round == mc_rounds){
            select_j <- T
          }
        }
      }

    }

    if(select_j){
      return(j)
    } else{
      return(NULL)
    }
  }

  selected <- union(selected, unlist(cali_selected))
  if(!is.null(X.pack$X.names))
    names(selected) <- X.pack$X.names[selected]


  # candidates_org <- candidates
  # candidates <- sapply(candidates, function(candidate){
  #   prefer <- kn_prefer(candidate, kn_stats_obs, alpha)
  #   if(prefer) return(candidate)
  #   else return(NA)
  # })
  # candidates <- candidates[!is.na(candidates)]
  # missed <- intersect(setdiff(candidates_org, candidates), selected)
  # write.table(missed, paste0("missed-", alpha, ".csv"), sep = ",", col.names = FALSE, append = T)


  # predict the sign of each beta
  sign_predict <- rep(0, p)
  sign_predict[kn_selected] <- (sign(matrix(y, nrow = 1) %*% (X.pack$X - X.pack$X_kn)))[kn_selected]
  sign_predict[setdiff(selected, kn_selected)] <- sign(matrix(y, nrow = 1) %*% X.pack$vj_mat[, setdiff(selected, kn_selected)])

  # record the working parameter for recursive exploration
  if(!is.null(record) && record$iteration > 0){
    checked_so_far <- union(union(record$checked_so_far, candidates), selected)

    if(length(checked_so_far) == p){
      message("Has tested all the hypotheses via calibration.")
    }
  } else{
    checked_so_far <- union(candidates, selected)
  }

  record <- list(iteration = ifelse(is.null(record), 1, record$iteration+1),
                 statistic = statistic,
                 alpha = alpha,
                 kappa = kappa,
                 n_cores = n_cores,
                 kn.selected = kn_selected,
                 checked_so_far = checked_so_far,
                 next_check_num = min(p - length(checked_so_far), max(5, length(candidates))),
                 X.pack = X.pack)

  # prepare the result
  result <- structure(list(call = match.call(),
                           X = X,
                           Xk = X.pack$X_kn,
                           y = y,
                           kn.statistic = kn_stats_obs,
                           selected = selected,
                           sign_predict = sign_predict,
                           record = record),
                      class = 'cknockoff.result')

  return(result)
}


Rhat_to_rej_j <- function(alpha, kappa, n, p, pval_rej_set, kn_rej_set,
                          cali_with_kn){
  df <- n-p

  if(!is.null(pval_rej_set$left)){
    pval_rej_lb <- pt(pval_rej_set$left, df = df, lower.tail = T)
    pval_rej_ub <- pt(pval_rej_set$right, df = df, lower.tail = T)
    pval_rej_mass <- sum(pval_rej_ub - pval_rej_lb)
  } else{
    return(Inf)
  }

  if(!is.null(kn_rej_set$left)){
    kn_rej_lb <- pt(kn_rej_set$left, df = df, lower.tail = T)
    kn_rej_ub <- pt(kn_rej_set$right, df = df, lower.tail = T)
    kn_rej_mass <- sum(kn_rej_ub - kn_rej_lb)
  } else{
    return(Inf)
  }


  if(cali_with_kn){
    if((min(kn_rej_lb) - pval_rej_ub[1] > 0.1)
       || (pval_rej_lb[2] - max(kn_rej_ub)) > 0.1){
      pval_rej_mass <- pval_rej_mass / 2
    }
  }

  b_j_est <- kappa * alpha / (ceiling(1/alpha - 1)) * (pval_rej_mass + kn_rej_mass) + (1-kappa) * alpha / p
  Rhat_to_rej <- pval_rej_mass / b_j_est

  return(Rhat_to_rej)
}

cknockoff_Rhat <- function(X.pack, y, j_exclude,
                           kn_stats_obs, tvals_obs,
                           statistic,
                           alpha,
                           kappa,
                           cali_with_kn,
                           Rhat_max_num,
                           Rhat_calc_max_step){
  n <- length(y)
  p <- length(tvals_obs)
  df <- n - p

  pvals_obs <- pvals_t(tvals_obs, df, side = "two")

  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha,
                               record = NULL, selected = NULL, cali_with_kn)
  candidates <- intersect(candidates, which(pvals_obs <= min(0.5 * alpha / p, 1e-3)))
  candidates <- setdiff(candidates, j_exclude)
  rank_p <- rank(pvals_obs[candidates])
  rank_kn <- rank(-abs(kn_stats_obs[candidates]))
  ranks <- rank(rank_p + rank_kn, ties.method = "first")
  candidates[ranks] <- candidates

  n_candidates <- min(length(candidates), Rhat_max_num)
  if(n_candidates > 0)  candidates <- candidates[1:n_candidates]
  else return(list(selected = NULL))

  y.pack <- process_y(X.pack, y)

  selected <- sapply(candidates, function(j){
    if(pvals_obs[j] < 1e-14) return(j)

    step_size <- pvals_obs[j] / (kappa * alpha / (ceiling(1/alpha - 1))) * 1.5
    pval_left_obs <- pt(tvals_obs[j], df = df, lower.tail = TRUE)

    calc_steps_remains <- Rhat_calc_max_step
    DP_j_accumulate <- pvals_obs[j]
    b_j_accumulate <- (1-kappa) * alpha / p

    pval_left_nodes <- pval_left_obs - sign(tvals_obs[j]) * (1:Rhat_calc_max_step) * step_size
    tval_nodes <- qt(pval_left_nodes, df = df, lower.tail = TRUE)
    vjy_nodes <- tj_to_vjy(tval_nodes, y.pack$y_Pi_Xnoj_res_norm2, df)
    y_results <- vjy_to_y(vjy_nodes, j, y.pack,
                          X.pack = list(vj = X.pack$vj_mat[, j],
                                        X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                        XXk_res_unit = X.pack$XXk_res_unit))

    for(mc_i in 1:Rhat_calc_max_step){
      if("sigma_tilde" %in% names(formals(statistic))){
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i],
                                sigma_tilde = y_results$sigmahat_XXk_res[mc_i])
      } else{
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i])
      }

      fj_result <- calc_fj(j, y_results$y_samples[, mc_i], sample_weight_mc = 1,
                           X.pack = NA, X.pack$vj_mat,
                           alpha, kappa, pvals_obs[j],
                           statistic, kn_stat_mc, y_results$sigmahat_X_res[mc_i],
                           cali_with_kn, Rhat_level = 0,
                           Rhat_max_num, Rhat_calc_max_step)

      DP_j_accumulate <- DP_j_accumulate + fj_result$DP_j * step_size
      b_j_accumulate <- b_j_accumulate + fj_result$b_j * step_size

      if(DP_j_accumulate <= b_j_accumulate) return(j)

      calc_steps_remains <- calc_steps_remains - 1
      if(calc_steps_remains == 0) break
    }

    if(calc_steps_remains > 0){
      pval_left_nodes <- 1-pval_left_obs + sign(tvals_obs[j]) * (1:Rhat_calc_max_step) * step_size
      tval_nodes <- qt(pval_left_nodes, df = df, lower.tail = TRUE)
      vjy_nodes <- tj_to_vjy(tval_nodes, y.pack$y_Pi_Xnoj_res_norm2, df)
      y_results <- vjy_to_y(vjy_nodes, j, y.pack,
                            X.pack = list(vj = X.pack$vj_mat[, j],
                                          X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                          XXk_res_unit = X.pack$XXk_res_unit))

      for(mc_i in 1:Rhat_calc_max_step){
        if("sigma_tilde" %in% names(formals(statistic))){
          kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i],
                                  sigma_tilde = y_results$sigmahat_XXk_res[mc_i])
        } else{
          kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i])
        }

        fj_result <- calc_fj(j, y_results$y_samples[, mc_i], sample_weight_mc = 1,
                             X.pack = NA, X.pack$vj_mat,
                             alpha, kappa, pvals_obs[j],
                             statistic, kn_stat_mc, y_results$sigmahat_X_res[mc_i],
                             cali_with_kn, Rhat_level = 0,
                             Rhat_max_num, Rhat_calc_max_step)


        DP_j_accumulate <- DP_j_accumulate + fj_result$DP_j * step_size
        b_j_accumulate <- b_j_accumulate + fj_result$b_j * step_size

        if(DP_j_accumulate <= b_j_accumulate) return(j)

        calc_steps_remains <- calc_steps_remains - 1
        if(calc_steps_remains == 0) break
      }

      return(NA)
    }
  })

  selected <- selected[!is.na(selected)]

  return(list(selected = selected))
}


calc_fj <- function(j, y_mc, sample_weight_mc, X.pack, vj_mat, alpha, kappa,
                    pval_obs, statistic, kn_stat_mc, sigmahat_X_res_mc,
                    cali_with_kn = F, Rhat_level = 0,
                    Rhat_max_num, Rhat_calc_max_step,
                    cali_var, cali_stat_obs, cali_stat_mc){
  n <- length(y_mc)
  p <- length(kn_stat_mc)
  df <- n - p

  ## compute weights
  # make rejection and fdp estimation based on knockoff statistics
  eskn_result <- kn.select(kn_stat_mc, kappa * alpha,
                         selective = T, early_stop = T)
  eskn_selected_mc <- eskn_result$selected

  # compute weight
  kn_null_num <- max(1, (eskn_result$fdp_est * length(eskn_selected_mc)))
  b_j <- kappa * alpha * (j %in% eskn_selected_mc) / kn_null_num

  ## compute individual FDP contribution (DP_j in the code)
  # p-values
  tvals_mc <- y_to_t(y_mc, vj_mat, sigmahat_X_res_mc)
  pvals_mc <- pvals_t(tvals_mc, df, side = "two")

  # Bonferroni rejections
  Bonf_selected_mc <- which(pvals_mc <= (1-kappa) * alpha / p)

  # knockoff rejections
  kn_result <- kn.select(kn_stat_mc, kappa * alpha,
                         selective = T, early_stop = F)
  kn_selected_mc <- kn_result$selected

  # R hat
  selected_mc <- union(kn_selected_mc, Bonf_selected_mc)

  # compute selections by a marginal p-value test
  # cali_selected_mc <- ifelse(pvals_mc[j] <= pval_obs, j, 0)
  # cali_selected_mc <- ifelse((pvals_mc[j] <= pval_obs) &&
  #                            kn_prefer(j, kn_stat_mc, alpha), j, 0)

  if(cali_var == "pval"){
    if(cali_with_kn){
      cali_selected_mc <- ifelse((pvals_mc[j] <= pval_obs) &&
                                   kn_prefer(j, kn_stat_mc, alpha), j, 0)
    } else{
      cali_selected_mc <- ifelse(pvals_mc[j] <= pval_obs, j, 0)
    }
  } else{
    cali_selected_mc <- ifelse(cali_stat_mc >= cali_stat_obs, j, 0)
  }

  if(Rhat_level > 0 && length(union(selected_mc, j)) == 1 &&
     (j %in% union(kn_selected_mc, cali_selected_mc))){
    Rhat_recursive <- cknockoff_Rhat(X.pack, y_mc, j,
                                     kn_stat_mc, tvals_mc,
                                     statistic,
                                     alpha,
                                     kappa,
                                     cali_with_kn = cali_with_kn,
                                     Rhat_max_num,
                                     Rhat_calc_max_step)$selected

    selected_mc <- union(selected_mc, Rhat_recursive)
    Rhat <- length(union(selected_mc, j))
  } else{
    Rhat <- NULL
  }

  DP_j <- (j %in% union(kn_selected_mc, cali_selected_mc)) / length(union(selected_mc, j))

  fj <- (DP_j - b_j) * sample_weight_mc - (1-kappa) * alpha / p

  return(list(fj = fj, DP_j = DP_j, b_j = b_j, Rhat = Rhat))
}


cali_statasitic <- function(X, X_kn, y, statistic, sigma_tilde, y_fit_noj_all,
                            type = c("W_abs", "beta_lasso", "LM_lasso")){
  n <- NROW(X)
  p <- NCOL(X)
  if(type == "W_abs"){
    cali_stat <- abs(statistic(X, X_kn, y, sigma_tilde))
  } else if(type == "beta_lasso"){
    lambda <- 2 * sigma_tilde / n
    fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=T, standardize=F, standardize.response=F, family="gaussian")
    cali_stat <- abs(as.vector(fit$beta))
  } else if(type == "LM_lasso"){
    cali_stat <- rep(NA, p)
    for(j in 1:p){
      y_res <- c(y - c(y_fit_noj_all[, j]))
      cali_stat[j] <- abs(sum(X[, j] * y_res))
    }
  }

  return(cali_stat)
}








