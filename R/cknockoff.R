#' @import stats methods
NULL


#' The cKnockoff procedure
#'
#' This function applies the cKnockoff procedure to the data \eqn{(X, y)}
#' subjecting to the Gaussian linear model, selecting variables/features
#' relevant for predicting the outcome with FDR control.
#'
#' @param X n-by-p matrix or data frame of features.
#' @param y response vector of length n.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero
#' (FALSE).
#' @param knockoffs method used to construct the knockoff matrix for X.
#' It should be a function taking a n-by-p matrix as input and returning a n-by-p
#' matrix of knockoff variables.
#' It is the same as the \code{knockoffs} argument in
#' \code{knockoff::knockoff.filter}.
#' By default, \code{ckn.create.fixed} is used.
#' @param statistic the knockoff feature statistics (W-statistics) function used
#' to assess variable importance.
#' Any function in the family "statistics" in the R package \code{knockoff} that
#' are suitable for the fixed-X setting can be supplied to this argument.
#' But please be aware of the efficiency issue as this function will be called
#' repeatedly in cknockoff. We suggest use the statistic functions provided by
#' our package.
#' By default, a lasso statistic from our package is used.
#' @param alpha  nominal false discovery rate (default: 0.05).
#' @param Rstar_refine A logical value determining if we use a better estimation
#' of the number of rejections in calibration.
#' If \code{TRUE}, the procedure is cKnockoff* and otherwise is the vanilla
#' cKnockoff.
#' The default is \code{FALSE}.
#' @param n_cores the number of cores to be used in computing cKnockoff in parallel.
#' package \code{doParallel} is required if \code{n_cores} > 1.
#' Otherwise it's computed sequentially.
#' @param prelim_result either a \code{knockoff.result} object returned by
#' \code{knockoff::knockoff.filter} or a \code{cknockoff.result} object returned
#' by \code{cknockoff}.
#' cknockoff can read the information from a knockoff result or a cknockoff
#' result and make possibly more rejections with the same FDR control.
#' See more details below.
#'
#' @param X.pack An object of class \code{cknockoff.X.pack} returned by function
#' \code{process_X}.
#' This is used only for simulations studies, where cknockoff is applied many
#' times to the same fixed X, to accelerate computation.
#' General users should ignore it and leave it as default = NULL.
#' If X.pack is provided, the argument "knockoffs" is not needed and will be
#' overwritten.
#'
#' @return An object of class \code{cknockoff.result}. It is a list similar
#' to the \code{knockoff.result} object, containing essentially the same
#' information:
#'  \item{X}{matrix of original variables (scaled and possibly augmented)}
#'  \item{Xk}{matrix of knockoff variables (cooresponding to the returned X)}
#'  \item{y}{response vector (possibly augmented)}
#'  \item{kn.statistic}{the knockoff feature statistics}
#'  \item{selected}{named vector of selected variables}
#'  \item{sign_predict}{the predicted signs of beta}
#'  \item{record}{a list recording some information during the calculation. It
#'  aims to make computing based on "prelim_result" possible and easy.
#'  Users may ignore it.}
#'
#' @details
#'
#' If argument \code{prelim_result} is supplied with a object, all the other
#' parameters except \code{n_cores} and \code{Rstar_refine} will be overwritten
#' by the information retrieved from it.
#'
#' If a \code{knockoff.result} object is supplied to \code{prelim_result},
#' cknockoff will return possibly more selections on the same problem with the
#' same FDR control.
#' To use, the function call of \code{knockoff::knockoff.filter} that returned
#' such \code{knockoff.result} object must have arguments "statistic"
#' and "fdr" explicitly provided (rather than relying on their defaults).
#' See examples below.
#'
#' It's possible that cknockoff cannot fetch the function
#' "statistic" or value for "alpha" by their names from a
#' \code{knockoff.result} object.
#' If this is the case, please supply them to cknockoff explicitly via arguments.
#'
#' It's possible to even make more rejections based on a \code{cknockoff.result}
#' object. This is because cknockoff will only explore those promising features
#' and discards the others, say those with large p-values.
#' Calling cknockoff() on a \code{cknockoff.result} object would further
#' explore some most promising ones among the discarded features and (rarely)
#' possibly make several more rejections.
#' Users can do this recursively and FDR is proved to be controlled whenever
#' they decide to stop.
#' The computational time of each call should be similar.
#'
#' If the \code{cknockoff.result} object is obtained by setting
#' \code{Rstar_refine = FALSE}, more rejections may be made with the same FDR
#' control by supplying \code{cknockoff.result} to \code{prelim_result} and set
#' \code{Rstar_refine = TRUE}. See examples below.
#'
#' @references
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 2.5 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#' print(which(1:p %in% nonzero))
#'
#' # Basic usage
#' result <- cknockoff(X, y, alpha = 0.05, n_cores = 1)
#' print(result$selected)
#'
#' # knockoff rejection
#' library("knockoff")
#' kn.result <- knockoff.filter(X, y,
#'                              knockoffs = ckn.create.fixed,
#'                              statistic = stat.glmnet_coefdiff_lm,
#'                              # must specify this argument explicitly
#'                              fdr = 0.05
#'                              # must specify this argument explicitly
#'                              )
#' print(kn.result$selected)
#'
#' # improve knockoff result
#' result <- cknockoff(prelim_result = kn.result, n_cores = 2)
#' print(result$selected)
#'
#' # improve previous cknockoff result
#' result <- cknockoff(prelim_result = result)
#' print(result$selected)
#'
#' # improve previous cknockoff result with cknockoff*
#' result <- cknockoff(prelim_result = result, n_cores = 2, Rstar_refine = TRUE)
#' print(result$selected)
#'
#' @export
cknockoff <- function(X, y,
                      intercept = TRUE,
                      knockoffs = ckn.create.fixed,
                      statistic = stat.glmnet_coefdiff_lm,
                      alpha = 0.05,
                      Rstar_refine = FALSE,
                      n_cores = 1,
                      prelim_result = NULL,
                      X.pack = NULL){

  mc_rounds <- 5
  mc_size <- 100  # must be even as paring samples is employed

  Rstar_max_try <- 3
  Rstar_calc_max_step <- 3

  # knockoff.type <- match.arg(knockoff.type)

  # cknockoff.call <- match.call.defaults()
  # cknockoff.call[[1]] <- as.symbol("parse_args")
  # args <- eval(cknockoff.call)
  envir <- parent.frame()
  args <- parse_args(X, y, intercept,
                     knockoffs, statistic,
                     alpha, Rstar_refine,
                     n_cores,
                     prelim_result, X.pack,
                     envir = envir)

  check_args(args)

  args <- process_args(args)
  X.pack <- args$X.pack
  y.pack <- args$y.pack
  intercept <- args$intercept
  X <- X.pack$X
  y <- y.pack$y.data$y
  statistic <- args$statistic
  alpha <- args$alpha
  n_cores <- args$n_cores
  Rstar_refine <- args$Rstar_refine
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

  # knockoff rejection
  if(is.null(record)){
    if("sigma_tilde" %in% names(formals(statistic))){
      kn_stats_obs <- statistic(X, X.pack$X_kn, y,
                                sigma_tilde = sqrt(y.pack$y.data$RSS_XXk
                                                   / y.pack$y.data$df_XXk))
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
  tvals_obs <- y_to_t(y, X.pack$vj_mat, sqrt(y.pack$y.data$RSS_X
                                             / y.pack$y.data$df_X))
  pvals_obs <- pvals_t(tvals_obs, df = y.pack$y.data$df_X, side = "two")

  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha, record, selected)

  # compute the observed calibration statistics
  for(j in candidates){
    sigma_est <- sqrt(y.pack$RSS_Xnoj[j] / (y.pack$y.data$df_X+1)) # overestimate sigma
    fit_on_rest <- glmnet::glmnet(X[, -j], y,
                                  lambda = 2 * sigma_est / n,
                                  intercept=T, standardize=F, standardize.response=F,
                                  family="gaussian")
    y.pack$Xy_bias[j] <- t(X[, j]) %*% (X[, -j] %*% fit_on_rest$beta + rep(fit_on_rest$a0, each = n))
  }
  cali_stats_obs <- abs(c(matrix(y, nrow = 1) %*% X) - y.pack$Xy_bias)

  # the confidence sequence alpha for controlling Monte-Carlo error
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
    sample_coupling <- couple_samples(y.pack$y.data$df_X, cali_rej_set, kn_rej_set)

    # storage for recording the calculated statistics for post-hoc Rstar refinement
    calc_Rstar_online <- F
    record_mc <- list(kn_stat = matrix(NA, nrow = p, ncol = mc_size),
                      DPj = rep(NA, mc_size), bj = rep(NA, mc_size))

    # divide the whole MC samples set into "mc_rounds" batches of size "mc_size"
    # and work on each batch sequentially for efficiency reason.
    for(mc_round in 1:mc_rounds){
      # expand the storage
      Ej_mc <- c(Ej_mc, rep(NA, mc_size))
      Ej_samples <- c(Ej_samples, rep(NA, mc_size))

      # generate Monte-Carlo samples of y conditional Sj for this batch
      samples.res <- y_sampler_cond_Sj(mc_size, cali_rej_set, kn_rej_set, j,
                                       y.pack, X.pack, sample_coupling)

      # make decision if we can tell whether Ej <=0 based on the sampling region
      if(is.numeric(samples.res) && mc_round == 1){
        if(samples.res == 0) break
        else if(samples.res == 1){
          break
        }
        else if(samples.res == -1){
          select_j <- T
          break
        }
      }

      # retrive the samples
      samples.y <- samples.res$samples.y
      samples.weight <- samples.res$samples.weight
      sigmahat.XXk_res <- sqrt(samples.res$samples.RSS_XXk / samples.res$df_XXk)

      # the upper and lower bound of Ej for constructing confidence sequence.
      Ej_bounds <- get_Ej_bound(alpha, p, samples.res$weights, sample_coupling)

      # for each y MC sample in conditional calibration
      for(mc_i in 1:mc_size){

        # compute knockoff stat
        if("sigma_tilde" %in% names(formals(statistic))){
          kn_stat_mc <- statistic(X, X.pack$X_kn, samples.y[, mc_i],
                                  sigma_tilde = sigmahat.XXk_res[mc_i])
        } else{
          kn_stat_mc <- statistic(X, X.pack$X_kn, samples.y[, mc_i])
        }

        # compute calibration stat
        cali_stat_mc <- abs(sum(X[, j] * samples.y[, mc_i]) - y.pack$Xy_bias[j])

        # do Rstar refinement in realtime?
        if(calc_Rstar_online){
          # prepare data for computing Rstar
          y.data_mc <- list(y = samples.y[, mc_i],
                            RSS_X = samples.res$samples.RSS_X[mc_i],
                            RSS_XXk = samples.res$samples.RSS_XXk[mc_i],
                            df_X = samples.res$df_X,
                            df_XXk = samples.res$df_XXk)
          Rstar.pack <- list(X.pack = X.pack, y.data = y.data_mc,
                            statistic = statistic,
                            Rstar_max_try = Rstar_max_try,
                            Rstar_calc_max_step = Rstar_calc_max_step)
        } else{
          Rstar.pack <- NA
        }
        # compute fj from the current MC sample
        fj_result <- calc_fj(j,
                             alpha, kn_stat_mc,
                             cali_stats_obs[j], cali_stat_mc,
                             Rstar.pack)

        # record the raw result
        mc_used <- mc_used + 1
        Ej_mc[mc_used] <- fj_result$fj * samples.weight[mc_i]

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

        # record the more detailed statistics for post-hoc Rstar refinement
        if(Rstar_refine && mc_round == 1){
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

      # post-hoc Rstar refinement, for those not rejected only
      if(Rstar_refine && !decision$reject && mc_round == 1){
        # prepare data for computing Rstar
        Rstar.pack <- list(X.pack = X.pack,
                           statistic = statistic,
                           Rstar_max_try = Rstar_max_try,
                           Rstar_calc_max_step = Rstar_calc_max_step)

        refined_result <- post_refine_Rstar(samples.res, record_mc,
                                            mc_used, n_Ej_samples,
                                            Ej_mc, Ej_samples,
                                            Ej_bounds, sample_coupling,
                                            decision, j,
                                            Rstar.pack,
                                            alpha, rej_alpha)

        decision <- refined_result$decision
        calc_Rstar_online <- refined_result$calc_Rstar_online
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
                 X.pack = X.pack,
                 y.pack = y.pack,
                 statistic = statistic,
                 alpha = alpha,
                 n_cores = n_cores,
                 kn.selected = kn_selected,
                 checked_so_far = checked_so_far,
                 next_check_num = min(p - length(checked_so_far), max(5, length(candidates))),
                 Rstar_refine = Rstar_refine)

  # prepare the result
  result <- structure(list(call = match.call(),
                           X = X.pack$X.org,
                           Xk = X.pack$X_kn.org,
                           y = y.pack$y.data$y.org,
                           intercept = intercept,
                           kn.statistic = kn_stats_obs,
                           selected = selected,
                           sign_predict = sign_predict,
                           record = record),
                      class = 'cknockoff.result')

  return(result)
}




## inner functions that compute important quantities in cknockoff()

# Replace some R^kn by R^* to do cKnockoff* after computing the f_j for cKnockoff,
# see Appendix of the cKnockoff paper for details
post_refine_Rstar <- function(samples.res, record_mc,
                              mc_used, n_Ej_samples,
                              Ej_mc, Ej_samples,
                              Ej_bounds, sample_coupling,
                              decision, j,
                              Rstar.pack,
                              alpha, rej_alpha){

  # make code shorter
  record_mc$DPj <- record_mc$DPj[1:mc_used]
  record_mc$bj <- record_mc$bj[1:mc_used]
  samples.weight <- samples.res$samples.weight[1:mc_used]
  # locate the samples with DP_j = 1, where Rstar refinement is needed
  Rstar_index <- which(abs(record_mc$DPj - 1) < 1e-10)

  if(length(Rstar_index) > 0 && length(Rstar_index) < mc_used){
    # how large should Rstar be to make j rejected by cKnockoff
    Rstar_to_rej_j <- calc_Rstar_to_rej_j(samples.weight, record_mc, Rstar_index)

    if(Rstar_to_rej_j <= Rstar.pack$Rstar_max_try){ # if it is within the budget
      Rstar_vec <- NULL # record the realized refined Rstar

      for(mc_i in Rstar_index){
        y.data_mc <- list(y = samples.res$samples.y[, mc_i],
                          RSS_X = samples.res$samples.RSS_X[mc_i],
                          RSS_XXk = samples.res$samples.RSS_XXk[mc_i],
                          df_X = samples.res$df_X,
                          df_XXk = samples.res$df_XXk)
        # compute Rstar
        Rstar <- cknockoff_Rstar(Rstar.pack$X.pack,
                                 y.data_mc,
                                 j_exclude = j,
                                 kn_stats_obs = record_mc$kn_stat[, mc_i],
                                 statistic = Rstar.pack$statistic,
                                 alpha = alpha,
                                 Rstar_max_try = Rstar.pack$Rstar_max_try,
                                 Rstar_calc_max_step = Rstar.pack$Rstar_calc_max_step)$Rstar

        # update the record for Ej
        Ej_mc[mc_i] <- (1/Rstar - record_mc$bj[mc_i]) * samples.weight[mc_i]
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

        # stop refining Rstar if we cannot make Rstar as large as Rstar_to_rej_j
        Rstar_vec <- c(Rstar_vec, Rstar)
        if(length(Rstar_vec) >= 10){
          if(mean(1/Rstar_vec) >= 1.5/Rstar_to_rej_j && !decision$reject){
            break
          }
        }

      }

      # after Rstar refinement is done, make decision based on all samples
      if(!decision$confident){
        decision <- make_decision(Ej_samples[1:n_Ej_samples], Ej_bounds, threshold = 0,
                                  rej_alpha = rej_alpha, accept_alpha = 0.05)
      }
    }
  }

  # do Rstar refinement online in the feture calculation,
  # make effect only when !decision$confident and decision$reject
  calc_Rstar_online <- T

  return(list(decision = decision, calc_Rstar_online = calc_Rstar_online,
              Ej_mc = Ej_mc, Ej_samples = Ej_samples))
}

# estimate the smallest number of R^* to make j rejected
calc_Rstar_to_rej_j <- function(samples.weight, record_mc, Rstar_index){
  # the values in fj that Rstar refinement cannot change
  other_value <- sum(record_mc$DPj[-Rstar_index] * samples.weight[-Rstar_index]) - sum(record_mc$bj * samples.weight)
  if(other_value >= 0){ # if fj > 0 even when Rstar = Inf
    Rstar_to_rej <- Inf
  } else{  # compute the Rstar needed to turn fj < 0
    Rstar_to_rej <- 1/(-other_value / sum(samples.weight[Rstar_index]))
  }

  return(Rstar_to_rej)
}

# compute R^* efficiently
cknockoff_Rstar <- function(X.pack, y.data, j_exclude,
                            kn_stats_obs,
                            statistic,
                            alpha,
                            Rstar_max_try,
                            Rstar_calc_max_step){
  y <- y.data$y

  n <- NROW(X.pack$X)
  p <- NCOL(X.pack$X)
  df_X <- y.data$df_X

  # p-values for screening
  tvals_obs <- y_to_t(y, X.pack$vj_mat, sqrt(y.data$RSS_X / df_X))
  pvals_obs <- pvals_t(tvals_obs, df_X, side = "two")
  # screening
  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha,
                               record = NULL, selected = NULL)
  # only keep the hypos with very small p-values
  candidates <- intersect(candidates,
                          which(pvals_obs <=
                                  min(alpha / p, 0.01 * (alpha / (ceiling(1/alpha - 1))) )))
  # exclude j as j must be in the Rstar
  candidates <- setdiff(candidates, j_exclude)
  # order the candidates by how promising they are
  rank_p <- rank(pvals_obs[candidates])
  rank_kn <- rank(-abs(kn_stats_obs[candidates]))
  ranks <- rank(rank_p + rank_kn, ties.method = "first")
  candidates[ranks] <- candidates

  # only keep the Rstar_max_try most promising hypos
  n_candidates <- min(length(candidates), Rstar_max_try)
  if(n_candidates > 0)  candidates <- candidates[1:n_candidates]
  else return(list(selected = NULL, Rstar = 1))

  # prepare data
  y.pack <- process_y(X.pack, y.data)

  for(j in candidates){
    sigma_est <- sqrt(y.pack$RSS_Xnoj[j] / (y.pack$y.data$df_X+1)) # overestimate sigma
    fit_on_rest <- glmnet::glmnet(X.pack$X[, -j], y,
                                  lambda = 2 * sigma_est / n,
                                  intercept=T, standardize=F, standardize.response=F,
                                  family="gaussian")
    y.pack$Xy_bias[j] <- t(X.pack$X[, j]) %*% (X.pack$X[, -j] %*% fit_on_rest$beta +
                                                 rep(fit_on_rest$a0, each = n))
  }
  cali_stats_obs <- abs(c(matrix(y, nrow = 1) %*% X.pack$X) - y.pack$Xy_bias)

  # decide if each hypo can be rejected
  selected <- sapply(candidates, function(j){
    # the calibration rejection set
    cali_rej_set <- where_cali_rej(j, y.pack, X.pack)

    if(is.null(cali_rej_set$left)) return(j)

    cali_rej_lb <- pt(cali_rej_set$left, df = df_X, lower.tail = T)
    cali_rej_ub <- pt(cali_rej_set$right, df = df_X, lower.tail = T)
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
    pval_left_start <- pt(tval_start, df = df_X, lower.tail = TRUE)

    # the DP_j and b_j value in the numerical integration up to now
    DP_j_accumulate <- cali_rej_mass
    b_j_accumulate <- 0

    # move grid point forward
    pval_left_nodes <- pval_left_start + direction * (1:Rstar_calc_max_step) * step_size
    # convert the p-value point to a y sample
    tval_nodes <- qt(pval_left_nodes, df = df_X, lower.tail = TRUE)
    vjy_nodes <- tj_to_vjy(tval_nodes, y.pack$RSS_Xnoj, df_X)

    y_results <- vjy_to_y(vjy_nodes, j, y.pack, X.pack)
    y_nodes <- y_results$y
    sigmahat.XXk_res <- sqrt(y_results$RSS_XXk / y_results$df_XXk)

    # at each grid point, compute fj
    for(mc_i in 1:Rstar_calc_max_step){
      if("sigma_tilde" %in% names(formals(statistic))){
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, mc_i],
                                sigma_tilde = sigmahat.XXk_res[mc_i])
      } else{
        kn_stat_mc <- statistic(X.pack$X, X.pack$X_kn, y_nodes[, mc_i])
      }

      cali_stat_mc <- abs(sum(X.pack$X[, j] * y_nodes[, mc_i]) - y.pack$Xy_bias[j])

      # compute fj without further Rstar refinement
      fj_result <- calc_fj(j,
                           alpha, kn_stat_mc,
                           cali_stats_obs[j], cali_stat_mc,
                           Rstar.pack = NA)

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

  return(list(selected = selected, Rstar = length(union(selected, j_exclude))))
}

# compute f_j
calc_fj <- function(j,
                    alpha, kn_stat_mc,
                    cali_stat_obs, cali_stat_mc,
                    Rstar.pack){

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
  # Rstar refinement?
  if(!is.na(Rstar.pack) && length(union(selected_mc, j)) == 1 &&
     (j %in% union(kn_selected_mc, cali_selected_mc))){
    Rstar_recursive <- cknockoff_Rstar(X.pack = Rstar.pack$X.pack,
                                       y.data = Rstar.pack$y.data,
                                       j_exclude = j,
                                       kn_stats_obs = kn_stat_mc,
                                       statistic = Rstar.pack$statistic,
                                       alpha = alpha,
                                       Rstar_max_try = Rstar.pack$Rstar_max_try,
                                       Rstar_calc_max_step = Rstar.pack$Rstar_calc_max_step)$selected

    selected_mc <- union(selected_mc, Rstar_recursive)
    Rstar <- length(union(selected_mc, j))
  } else{
    Rstar <- NULL
  }

  # compute DP_j
  DP_j <- (j %in% union(kn_selected_mc, cali_selected_mc)) / length(union(selected_mc, j))

  # compute fj
  fj <- DP_j - b_j

  return(list(fj = fj, DP_j = DP_j, b_j = b_j, Rstar = Rstar))
}










