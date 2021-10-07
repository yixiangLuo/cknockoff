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
#' y = X %*% beta + rnorm(n)
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
                      X.pack = NULL){

  mc_rounds <- 5
  mc_size <- 100
  mc_size <- 2 * ceiling(mc_size / 2)

  # knockoff.type <- match.arg(knockoff.type)

  cknockoff.call <- match.call.defaults()
  cknockoff.call[[1]] <- as.symbol("parse_args")
  args <- eval(cknockoff.call)

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

    kn_selected <- kn_rej_fdp(kn_stats_obs, kappa * alpha,
                              selective = T, early_stop = 0)$selected
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
  candidates <- cKn_candidates(kn_stats_obs, pvals_obs, alpha, record, selected)

  # check each hypothesis in sequence/parallel
  cali_selected <- forall(j = candidates) %exec% {

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
    pval_rej_set <- where_pval_rej(j, y.pack)

    # we divide the whole MC samples set into "mc_rounds" batches of size "mc_size"
    # and work on each batch sequentially for efficiency reason.
    for(mc_round in 1:mc_rounds){
      ineq <- c(ineq, rep(NA, mc_size))
      ineq_combine <- c(ineq_combine, rep(NA, mc_size/2))

      # generate Monte-Carlo samples of y conditional Sj for this batch
      sample_res <- y_sampler_cond_Sj(mc_size, pval_rej_set, kn_rej_set, j,
                                      y.pack,
                                      X.pack = list(vj = X.pack$vj_mat[, j],
                                                    X_res_Xk_basis = X.pack$X_res_Xk_basis,
                                                    XXk_res_unit = X.pack$XXk_res_unit))
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

        ## compute weights
        # make rejection and fdp estimation based on knockoff statistics
        kn_result <- kn_rej_fdp(kn_stat_mc, kappa * alpha,
                                selective = T, early_stop = 1)
        kn_selected_RB <- kn_result$selected

        # compute weight
        kn_null_num <- max(1, (kn_result$fdp_est * length(kn_selected_RB)))
        w_tilde <- (j %in% kn_selected_RB) / kn_null_num * kappa

        ## compute individual FDP contribution (g_star in the code)
        # p-values
        tvals_RB <- y_to_t(y_cond[, mc_i], X.pack$vj_mat, sigmahat_X_res[mc_i])
        pvals_RB <- pvals_t(tvals_RB, df, side = "two")

        # Bonferroni rejections
        Bonf_selected_RB <- which(pvals_RB <= (1-kappa) * alpha / p)

        # knockoff rejections
        kn_result <- kn_rej_fdp(kn_stat_mc, kappa * alpha,
                                selective = T, early_stop = 0)
        kn_selected_RB <- kn_result$selected

        # R hat
        selected_RB <- union(kn_selected_RB, Bonf_selected_RB)

        # compute selections by a marginal p-value test
        Mg_selected_RB <- ifelse(pvals_RB[j] <= pval_obs, j, 0)

        g_star <- (j %in% union(kn_selected_RB, Mg_selected_RB)) / length(union(selected_RB, j))

        # record the result
        mc_used <- mc_used + 1
        ineq[mc_used] <- (g_star - w_tilde * alpha) * sample_weights[mc_i] - (1-kappa) * alpha / p

        # compute the confidence interval
        if(mc_used %% 2 == 0){

          combine_used <- mc_used/2
          ineq_combine[combine_used] <- mean(c(ineq[mc_used-1], ineq[mc_used]))

          # make decision
          decision <- make_decision(ineq_combine[1:combine_used], ineq_bounds, threshold = 0,
                                    rej_alpha = min(0.05, alpha * max(1, length(init_selected)) / length(candidates)),
                                    accept_alpha = 0.05)

          if(decision$confident){
            break
          }
        }
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

  # # predict the sign of each beta
  # sign_predict <- rep(0, p)
  # sign_predict[kn_selected] <- (sign(t(X - X.pack$X_kn) %*% y))[kn_selected]
  # sign_predict[setdiff(selected, kn_selected)] <- sign(matrix(y, nrow = 1) %*% X.pack$vj_mat[, setdiff(selected, kn_selected)])

  # record the working parameter for recursive exploration
  if(!is.null(record) && record$iteration > 0){
    checked_so_far <- union(union(record$checked_so_far, candidates), selected)
  } else{
    checked_so_far <- union(candidates, selected)
  }
  if(length(checked_so_far) == p){
    message("Has tested all the hypotheses via calibration.")
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
                           record = record),
                      class = 'cknockoff.result')

  return(result)
}
