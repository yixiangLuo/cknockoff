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
                      n_cores = 1,
                      knockoff.type = c("fixed", "model"),
                      prelim_result = NULL,
                      X.pack = NULL,
                      Rhat_refine = F){

  mc_rounds <- 5
  mc_size <- 100  # must be even as paring samples is employed

  Rhat_max_try <- 3
  Rhat_calc_max_step <- 3

  # knockoff.type <- match.arg(knockoff.type)

  # cknockoff.call <- match.call.defaults()
  # cknockoff.call[[1]] <- as.symbol("parse_args")
  # args <- eval(cknockoff.call)
  envir <- parent.frame()
  args <- parse_args(X, y, knockoffs, statistic,
                     alpha, n_cores, knockoff.type,
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

  # knockoff rejection
  if(is.null(record)){
    if("sigma_tilde" %in% names(formals(statistic))){
      kn_stats_obs <- statistic(X, X.pack$X_kn, y,
                                sigma_tilde = y.pack$sigmahat_XXk_res)
    } else{
      kn_stats_obs <- statistic(X, X.pack$X_kn, y)
    }

    kn_selected <- kn.select(kn_stats_obs, alpha,
                             selective = T, early_stop = F)$selected
  } else{
    kn_stats_obs <- record$kn.statistic
    kn_selected <- record$kn.selected
  }

  init_selected <- kn_selected

  if(!is.null(record) && record$iteration > 0){
    selected <- record$selected_so_far
  } else{
    selected <- init_selected
  }

  # find the promising features not selected or checked yet
  tvals_obs <- y_to_t(y, X.pack$vj_mat, sqrt(y.pack$y_Pi_X_res_norm2 / df))
  pvals_obs <- pvals_t(tvals_obs, df = df, side = "two")

  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha, record, selected)

  # compute the observed calibration statistics
  for(j in candidates){
    fit_on_rest <- glmnet::glmnet(X[, -j], y, lambda = 2 * y.pack$sigmahat_XXk_res/n, intercept=T, standardize=F, standardize.response=F, family="gaussian")
    y.pack$Xy_bias[j] <- t(X[, j]) %*% (X[, -j] %*% fit_on_rest$beta + rep(fit_on_rest$a0, each = n))
  }
  cali_stats_obs <- abs(c(matrix(y, nrow = 1) %*% X) - y.pack$Xy_bias)

  # the confidence sequence alpha for controling Monte-Carlo error
  rej_alpha <- min(0.05, alpha * max(1, length(init_selected)) / length(candidates))

  # check each hypothesis in sequence/parallel
  cali_selected <- forall(j = candidates, .options.multicore = list(preschedule = F)) %exec% {

    # prepare repeatedly used variables
    Ej_mc <- NULL  # raw record for the Ej samples
    mc_used <- 0  # number of MC samples used up to now
    Ej_samples <- NULL  # adjusted record for the Ej samples based on sample_coupling
    n_Ej_samples <- 0  # number of recorded Ej_samples up to now

    select_j <- F  # do we select j, initially set as F

    # prepare the sets where knockoff or the calibration stat would reject j
    # conditional on Sj. This will be used in generating Monte-Carlo samples of y
    # conditional on Sj.
    kn_rej_set <- where_kn_rej(alpha, j, y.pack, X.pack,
                               statistic, method = kn_region_method)

    cali_rej_set <- where_cali_rej(j, y.pack, X.pack)

    # couple the samples?
    sample_coupling <- couple_samples(n, p, cali_rej_set, kn_rej_set)

    # storage for recording the calculated statistics for post-hoc Rhat refinement
    calc_Rhat_online <- F
    record_mc <- list(kn_stat = matrix(NA, nrow = p, ncol = mc_size),
                      DPj = rep(NA, mc_size), bj = rep(NA, mc_size))

    # divide the whole MC samples set into "mc_rounds" batches of size "mc_size"
    # and work on each batch sequentially for efficiency reason.
    for(mc_round in 1:mc_rounds){
      # expand the storage
      Ej_mc <- c(Ej_mc, rep(NA, mc_size))
      Ej_samples <- c(Ej_samples, rep(NA, mc_size))

      # generate Monte-Carlo samples of y conditional Sj for this batch
      sample_res <- y_sampler_cond_Sj(mc_size, cali_rej_set, kn_rej_set, j,
                                      y.pack, X.pack, sample_coupling)

      # make decision if we can tell whether Ej <=0 based on the sampling region
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

      # the upper and lower bound of Ej for constructing confidence sequence.
      Ej_bounds <- get_Ej_bound(alpha, p, sample_res$weights, sample_coupling)

      # for each y MC sample in conditional calibration
      for(mc_i in 1:mc_size){

        # compute knockoff stat
        if("sigma_tilde" %in% names(formals(statistic))){
          kn_stat_mc <- statistic(X, X.pack$X_kn, y_cond[, mc_i],
                                  sigma_tilde = sigmahat_XXk_res[mc_i])
        } else{
          kn_stat_mc <- statistic(X, X.pack$X_kn, y_cond[, mc_i])
        }

        # compute calibration stat
        cali_stat_mc <- abs(sum(X[, j] * y_cond[, mc_i]) - y.pack$Xy_bias[j])

        # do Rhat refinement in realtime?
        if(calc_Rhat_online){
          Rhat.pack <- list(X.pack = X.pack, y_mc = y_cond[, mc_i],
                            statistic = statistic,
                            sigmahat_XXk_res = sigmahat_XXk_res[mc_i],
                            Rhat_max_try = Rhat_max_try,
                            Rhat_calc_max_step = Rhat_calc_max_step)
        } else{
          Rhat.pack <- NA
        }
        # compute fj from the current MC sample
        fj_result <- calc_fj(j,
                             alpha, kn_stat_mc,
                             cali_stats_obs[j], cali_stat_mc,
                             Rhat.pack)

        # record the raw result
        mc_used <- mc_used + 1
        Ej_mc[mc_used] <- fj_result$fj * sample_weights[mc_i]

        # record the adjusted Ej_samples based on coupling
        if(!sample_coupling){
          n_Ej_samples <- mc_used
          Ej_samples[n_Ej_samples] <- Ej_mc[mc_used]
          try_decide <- T
        } else{
          if(mc_used %% 2 == 0){
            n_Ej_samples <- mc_used/2
            Ej_samples[n_Ej_samples] <- (Ej_mc[mc_used-1] + Ej_mc[mc_used]) / 2
            try_decide <- T
          } else{
            try_decide <- F
          }
        }

        # record the more detailed statistics for post-hoc Rhat refinement
        if(Rhat_refine && mc_round == 1){
          record_mc$DPj[mc_i] <- fj_result$DP_j
          record_mc$bj[mc_i] <- fj_result$b_j
          if(abs(fj_result$DP_j - 1) < 1e-10){
            record_mc$kn_stat[, mc_i] <- kn_stat_mc
          }
        }

        # calculate the confidence interval and make decision
        if(try_decide){
          decision <- make_decision(Ej_samples[1:n_Ej_samples], Ej_bounds, threshold = 0,
                                    rej_alpha = rej_alpha, accept_alpha = 0.05)

          if(decision$confident){ break }
        }

      }

      # post-hoc Rhat refinement, for those not rejected only
      if(Rhat_refine && !decision$reject && mc_round == 1){
        Rhat.pack <- list(X.pack = X.pack,
                          statistic = statistic,
                          Rhat_max_try = Rhat_max_try,
                          Rhat_calc_max_step = Rhat_calc_max_step)

        refined_result <- post_refine_Rhat(sample_res, record_mc,
                                           mc_used, n_Ej_samples,
                                           Ej_mc, Ej_samples,
                                           Ej_bounds, sample_coupling,
                                           decision, j,
                                           Rhat.pack,
                                           alpha, rej_alpha)

        decision <- refined_result$decision
        calc_Rhat_online <- refined_result$calc_Rhat_online
        Ej_mc <- refined_result$Ej_mc
        Ej_samples <- refined_result$Ej_samples
      }

      # decide rejecting or not
      if(decision$confident){ # with confidence, do as suggested
        if(decision$reject){
          select_j <- T
        }
        break
      } else{
        # if suggest accepting without confidence, follow it (conservative)
        if(!decision$reject){
          break
        } else{ # if suggest rejecting without confidence, continue running more rounds
          if(mc_round == mc_rounds){ # if used up our budget, reject even without confidence
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

  # resemble the selection set
  selected <- union(selected, unlist(cali_selected))
  if(!is.null(X.pack$X.names))
    names(selected) <- X.pack$X.names[selected]

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




couple_samples <- function(n, p, cali_rej_set, kn_rej_set){
  df <- n-p

  # compute mass of the calibration rejection set
  rej_set <- cali_rej_set
  if(!is.null(rej_set$left)){
    cali_rej_lb <- pt(rej_set$left, df = df, lower.tail = T)
    cali_rej_ub <- pt(rej_set$right, df = df, lower.tail = T)
    rej_mass <- sum(cali_rej_ub - cali_rej_lb)
  } else{
    return(F)
  }

  # compute mass of the rest rejection set
  rest_set <- interval_minus(kn_rej_set, cali_rej_set)
  if(!is.null(rest_set$left)){
    cali_rest_lb <- pt(rest_set$left, df = df, lower.tail = T)
    cali_rest_ub <- pt(rest_set$right, df = df, lower.tail = T)
    rest_mass <- sum(cali_rest_ub - cali_rest_lb)
  } else{
    return(F)
  }

  # if rej_mass and rest_mass are not far from each other, couple the samples to reduce vairance
  # otherwise don't to avoid over sampling in a relatively small region
  threshold <- 3
  sample_coupling <- (rej_mass / rest_mass <= threshold) && (rest_mass / rej_mass <= threshold)

  return(sample_coupling)
}

post_refine_Rhat <- function(sample_res, record_mc,
                             mc_used, n_Ej_samples,
                             Ej_mc, Ej_samples,
                             Ej_bounds, sample_coupling,
                             decision, j,
                             Rhat.pack,
                             alpha, rej_alpha){

  # make code shorter
  record_mc$DPj <- record_mc$DPj[1:mc_used]
  record_mc$bj <- record_mc$bj[1:mc_used]
  sample_weights <- sample_res$sample_weights[1:mc_used]
  # locate the samples with DP_j = 1, where Rhat refinement is needed
  Rhat_index <- which(abs(record_mc$DPj - 1) < 1e-10)

  if(length(Rhat_index) > 0 && length(Rhat_index) < mc_used){
    # how large should Rhat be to make j rejected by cKnockoff
    Rhat_to_rej_j <- calc_Rhat_to_rej_j(sample_weights, record_mc, Rhat_index)

    if(Rhat_to_rej_j <= Rhat.pack$Rhat_max_try){ # if it is within the budget
      Rhat_vec <- NULL # record the realized refined Rhat

      for(mc_i in Rhat_index){
        # compute Rhat
        Rhat <- cknockoff_Rhat(Rhat.pack$X.pack,
                               sample_res$y_samples[, mc_i],
                               j_exclude = j,
                               kn_stats_obs = record_mc$kn_stat[, mc_i],
                               sigmahat_XXk_res = sample_res$sigmahat_XXk_res[mc_i],
                               statistic = Rhat.pack$statistic,
                               alpha = alpha,
                               Rhat_max_try = Rhat.pack$Rhat_max_try,
                               Rhat_calc_max_step = Rhat.pack$Rhat_calc_max_step)$Rhat

        # update the record for Ej
        Ej_mc[mc_i] <- (1/Rhat - record_mc$bj[mc_i]) * sample_weights[mc_i]
        # update the record for Ej_samples
        if(!sample_coupling){
          i_Ej_samples <- mc_i
          Ej_samples[i_Ej_samples] <- Ej_mc[mc_i]
        } else{
          i_Ej_samples <- ceiling(mc_i/2)
          Ej_samples[i_Ej_samples] <- (Ej_mc[2*i_Ej_samples-1] + Ej_mc[2*i_Ej_samples]) / 2
        }

        # make decision as if we have samples 1:mc_i
        decision <- make_decision(Ej_samples[1:i_Ej_samples], Ej_bounds, threshold = 0,
                                  rej_alpha = rej_alpha, accept_alpha = 0.05)
        if(decision$confident){ break }

        # stop refining Rhat if we cannot make Rhat as large as Rhat_to_rej_j
        Rhat_vec <- c(Rhat_vec, Rhat)
        if(length(Rhat_vec) >= 10){
          if(mean(1/Rhat_vec) >= 1.5/Rhat_to_rej_j && !decision$reject){
            break
          }
        }

      }

      # after Rhat refinement is done, make decision based on all samples
      if(!decision$confident){
        decision <- make_decision(Ej_samples[1:n_Ej_samples], Ej_bounds, threshold = 0,
                                  rej_alpha = rej_alpha, accept_alpha = 0.05)
      }
    }
  }

  # do Rhat refinement online in the feture calculation,
  # make effect only when !decision$confident and decision$reject
  calc_Rhat_online <- T

  return(list(decision = decision, calc_Rhat_online = calc_Rhat_online,
              Ej_mc = Ej_mc, Ej_samples = Ej_samples))
}

calc_Rhat_to_rej_j <- function(sample_weights, record_mc, Rhat_index){
  # the values in fj that Rhat refinement cannot change
  other_value <- sum(record_mc$DPj[-Rhat_index] * sample_weights[-Rhat_index]) - sum(record_mc$bj * sample_weights)
  if(other_value >= 0){ # if fj > 0 even when Rhat = Inf
    Rhat_to_rej <- Inf
  } else{  # compute the Rhat needed to turn fj < 0
    Rhat_to_rej <- 1/(-other_value / sum(sample_weights[Rhat_index]))
  }

  return(Rhat_to_rej)
}

cknockoff_Rhat <- function(X.pack, y, j_exclude,
                           kn_stats_obs,
                           sigmahat_XXk_res,
                           statistic,
                           alpha,
                           Rhat_max_try,
                           Rhat_calc_max_step){
  n <- length(y)
  p <- length(kn_stats_obs)
  df <- n - p

  # p-values for screening
  tvals_obs <- y_to_t(y, X.pack$vj_mat, sigmahat_XXk_res)
  pvals_obs <- pvals_t(tvals_obs, df, side = "two")
  # screening
  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha,
                               record = NULL, selected = NULL)
  # only keep the hypos with very small p-values
  candidates <- intersect(candidates, which(pvals_obs <= min(alpha / p,
                                                             0.01 * (alpha / (ceiling(1/alpha - 1))) )))
  # exclude j as j must be in the Rhat
  candidates <- setdiff(candidates, j_exclude)
  # order the candidates by how promising they are
  rank_p <- rank(pvals_obs[candidates])
  rank_kn <- rank(-abs(kn_stats_obs[candidates]))
  ranks <- rank(rank_p + rank_kn, ties.method = "first")
  candidates[ranks] <- candidates

  # only keep the Rhat_max_try most promising hypos
  n_candidates <- min(length(candidates), Rhat_max_try)
  if(n_candidates > 0)  candidates <- candidates[1:n_candidates]
  else return(list(selected = NULL, Rhat = 1))

  # prepare data
  y.pack <- process_y(X.pack, y)

  for(j in candidates){
    fit_on_rest <- glmnet::glmnet(X.pack$X[, -j], y, lambda = 2 * y.pack$sigmahat_XXk_res/n, intercept=T, standardize=F, standardize.response=F, family="gaussian")
    y.pack$Xy_bias[j] <- t(X.pack$X[, j]) %*% (X.pack$X[, -j] %*% fit_on_rest$beta + rep(fit_on_rest$a0, each = n))
  }
  cali_stats_obs <- abs(c(matrix(y, nrow = 1) %*% X.pack$X) - y.pack$Xy_bias)

  # decide if each hypo can be rejected
  selected <- sapply(candidates, function(j){
    # the calibration rejection set
    cali_rej_set <- where_cali_rej(j, y.pack, X.pack)

    if(is.null(cali_rej_set$left)) return(j)

    cali_rej_lb <- pt(cali_rej_set$left, df = df, lower.tail = T)
    cali_rej_ub <- pt(cali_rej_set$right, df = df, lower.tail = T)
    cali_rej_mass <- sum(cali_rej_ub - cali_rej_lb)

    # if cali_rej_set too large, don't reject j
    if(cali_rej_mass > 0.0067 * (alpha / (ceiling(1/alpha - 1)))) return(NA)

    # where to start the numerical integration grid point?
    if(length(cali_rej_set$left) == 1){
      # if cali_rej_set is an interval, start at the value closer to 0
      t_rej_points <- c(cali_rej_set$left, cali_rej_set$right)
      tval_start <- t_rej_points[which.min(abs(t_rej_points))]
    } else{
      # otherwise start at the observed value
      tval_start <- tvals_obs[j]
    }
    # forward the numerical integration towards 0
    direction <- -sign(tval_start)

    # step size in the numerical integration, which expects to reject j at one try
    step_size <- cali_rej_mass / (alpha / (ceiling(1/alpha - 1))) * 1.5
    # convert the t stat value to a one-sided p-value (so Lebesgue measure)
    pval_left_start <- pt(tval_start, df = df, lower.tail = TRUE)

    # the DP_j and b_j value in the numerical integration up to now
    DP_j_accumulate <- cali_rej_mass
    b_j_accumulate <- 0

    # move grid point forward
    pval_left_nodes <- pval_left_start + direction * (1:Rhat_calc_max_step) * step_size
    # convert the p-value point to a y sample
    tval_nodes <- qt(pval_left_nodes, df = df, lower.tail = TRUE)
    vjy_nodes <- tj_to_vjy(tval_nodes, y.pack$y_Pi_Xnoj_res_norm2, df)
    y_results <- vjy_to_y(vjy_nodes, j, y.pack, X.pack)

    # at each grid point, compute fj
    for(mc_i in 1:Rhat_calc_max_step){
      if("sigma_tilde" %in% names(formals(statistic))){
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i],
                                sigma_tilde = y_results$sigmahat_XXk_res[mc_i])
      } else{
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_results$y_samples[, mc_i])
      }

      cali_stat_mc <- abs(sum(X.pack$X[, j] * y_results$y_samples[, mc_i]) - y.pack$Xy_bias[j])

      # compute fj without further Rhat refinement
      fj_result <- calc_fj(j,
                           alpha, kn_stat_mc,
                           cali_stats_obs[j], cali_stat_mc,
                           Rhat.pack = NA)

      # update DP_j and b_j
      DP_j_accumulate <- DP_j_accumulate + fj_result$DP_j * step_size
      b_j_accumulate <- b_j_accumulate + fj_result$b_j * step_size

      # reject if we see Ej <= 0, as Ej is decreasing when grid point move forward
      if(DP_j_accumulate <= b_j_accumulate) return(j)
    }

    # if don't Ej <= 0, accept j
    return(NA)
  })

  # remove the NAs in the selection set
  selected <- selected[!is.na(selected)]

  return(list(selected = selected, Rhat = length(union(selected, j_exclude))))
}


calc_fj <- function(j,
                    alpha, kn_stat_mc,
                    cali_stat_obs, cali_stat_mc,
                    Rhat.pack){

  ## compute b_j term
  # make rejection and fdp estimation based on knockoff statistics
  eskn_result <- kn.select(kn_stat_mc, alpha,
                         selective = T, early_stop = T)
  eskn_selected_mc <- eskn_result$selected

  # compute b_j
  kn_null_num <- max(1, (eskn_result$fdp_est * length(eskn_selected_mc)))
  b_j <- alpha * (j %in% eskn_selected_mc) / kn_null_num

  ## compute DP_j term
  # knockoff rejections
  kn_result <- kn.select(kn_stat_mc, alpha,
                         selective = T, early_stop = F)
  kn_selected_mc <- kn_result$selected

  # calibration rejections
  cali_selected_mc <- ifelse(cali_stat_mc >= cali_stat_obs, j, 0)

  # R hat
  selected_mc <- kn_selected_mc
  # Rhat refinement?
  if(!is.na(Rhat.pack) && length(union(selected_mc, j)) == 1 &&
     (j %in% union(kn_selected_mc, cali_selected_mc))){
    Rhat_recursive <- cknockoff_Rhat(X.pack = Rhat.pack$X.pack,
                                     y = Rhat.pack$y_mc,
                                     j_exclude = j,
                                     kn_stats_obs = kn_stat_mc,
                                     sigmahat_XXk_res = Rhat.pack$sigmahat_XXk_res,
                                     statistic = Rhat.pack$statistic,
                                     alpha = alpha,
                                     Rhat_max_try = Rhat.pack$Rhat_max_try,
                                     Rhat_calc_max_step = Rhat.pack$Rhat_calc_max_step)$selected

    selected_mc <- union(selected_mc, Rhat_recursive)
    Rhat <- length(union(selected_mc, j))
  } else{
    Rhat <- NULL
  }

  # compute DP_j
  DP_j <- (j %in% union(kn_selected_mc, cali_selected_mc)) / length(union(selected_mc, j))

  # compute fj
  fj <- DP_j - b_j

  return(list(fj = fj, DP_j = DP_j, b_j = b_j, Rhat = Rhat))
}










