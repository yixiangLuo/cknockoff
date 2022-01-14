#' Creates fixed-X knockoff variables by MRC
#'
#' This function is a R wrapper of the \code{knocproduce_FX_knockoffskpy}
#' function in the Python package \code{knockpy}. It creates fixed-X knockoff
#' variables by Minimizing Reconstructability
#' (see \href{Spector and Janson (2021)}{https://arxiv.org/abs/2011.14625}).
#' To use this function, Users should have the R package \code{reticulate} and
#' the python package \code{knockpy} installed.
#'
#' @param X n-by-p matrix of the features.
#' @param method Same argument as the "method" in the \code{smatrix.compute_smatrix}
#' function in knockpy package.
#' Take value from "mvr", "maxent", "mmi", "sdp", "equicorrelated", and "ci".
#' The default is "mvr".
#' @param solver Same argument as the "solver" in the \code{smatrix.compute_smatrix}
#' function in knockpy package.
#' Take value from "cd", "sdp", "psgd".
#' Can be leave as \code{NULL}.
#' @param knockpy The R object wrapping the Python module "knockpy". Users may
#' provide their own wrapper using "reticulate::import" if "knockpy" is not
#' installed in the default python environment
#' (e.g. in some environments created by anaconda).
#' @param num_processes Same argument as the "num_processes" in the
#' \code{smatrix.compute_smatrix} function in knockpy package.
#' Number of parallel process to use in computing the matrix approximately.
#'
#' @return A list containing the following components:
#'  \item{X}{n-by-p matrix of original variables (normalized).}
#'  \item{Xk}{n-by-p matrix of knockoff variables.}
#'
#' @references
#'
#'
#' @export
create.fixed.MRC <- function(X,
                             method = c("mvr", "maxent", "mmi",
                                        "sdp", "equicorrelated", "ci"),
                             solver = NULL,
                             knockpy = NULL,
                             num_processes = 1){

  if(!requireNamespace("reticulate", quietly=T)) {
    stop("R package 'reticulate' is required for calling python package 'knockpy' but is not installed.")
  }

  method <- match.arg(method)
  if(is.null(solver)){
    solver <- ifelse(method == "sdp", "sdp", "cd")
  }

  if(is.null(knockpy)){
    tryCatch({
      knockpy <- reticulate::import("knockpy")
    }, error = function(msg){
      stop(msg)
    })
  }

  X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))

  Simga_X <- t(X) %*% X
  n <- NROW(X)

  s_mat <- knockpy$smatrix$compute_smatrix(Simga_X, method = method, solver = solver,
                                           num_processes = num_processes)
  Xk <- knockpy$knockoffs$produce_FX_knockoffs(X, solve(Simga_X), s_mat)
  Xk <- matrix(Xk, nrow = n)
  # Xk <- scale(Xk, center = FALSE, scale = sqrt(colSums(Xk^2)))

  return(list(X = X, Xk = Xk))
}
# @examples
# p <- 100; n <- 300; k <- 15
# X <- matrix(rnorm(n*p), n)
# nonzero <- sample(p, k)
# beta <- 2.5 * (1:p %in% nonzero)
# y <- X %*% beta + rnorm(n)
# print(which(1:p %in% nonzero))
#
# result <- cknockoff(X, y,
#                     knockoffs = create.fixed.MRC,
#                     alpha = 0.05, n_cores = 1)
# print(result$selected)


#' Select variables by knockoff.
#'
#' Select variables relevant for predicting the outcome of interest, using the
#' knockoff feature statistics and applying selective SeqStep.
#'
#' @param kn_statistics The knockoff feature statistics (W-statistics)
#' @param alpha target false discovery rate
#' @param selective If TRUE, use selective SeqStep for selection
#' (only select variables with a positive W-statistics);
#' if FALSE, use SeqStep for selection
#' (select variables if their W-statistics have large absolute values.)
#' @param early_stop If FALSE, the selection set is the same as in knockoff
#' (stop at \eqn{\tau});
#' If TRUE, the selection set is what we used to construct the budget \eqn{b_j}
#' (stop at \eqn{\tau_1}).
#'
#' @return A list containing the following components:
#'  \item{selected}{vector of selected variables.}
#'  \item{fdp_est}{the estimated FDP.}
#'  \item{W_k_hat}{The absolute value of the W-statistics where we stop.}
#'
#' @examples
#'
#' @keywords internal
kn.select <- function(kn_statistics, alpha,
                      selective = T, early_stop = F){
  p <- length(kn_statistics)
  W_desc_ind <- order(abs(kn_statistics), decreasing = T)
  W_desc <- kn_statistics[W_desc_ind]

  # fdp <- sapply(abs(W_desc), function(t){
  #     (1 + sum(kn_statistics <= -t)) / max(1, sum(kn_statistics >= t))
  # })
  fdp <- rep(NA, p)
  neg_stat_num <- 0
  pos_stat_num <- 0
  for(k in 1:p){
    pos_stat_num <- pos_stat_num + (W_desc[k] > 0)
    neg_stat_num <- neg_stat_num + (W_desc[k] < 0)
    fdp[k] <- (1 + neg_stat_num) / max(1, pos_stat_num)
  }
  ok <- which(fdp <= alpha)
  k_hat <- ifelse(length(ok) > 0, max(ok), 0)

  if(early_stop){
    ok <- which(cumsum(W_desc > 0) < 1/alpha)
    k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
  }
  # if(early_stop == 2){
  #   ok <- which((fdp <= 1) & (cumsum(W_desc > 0) < 1/alpha))
  #   k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
  # }

  if(k_hat > 0){
    selected <- sort(W_desc_ind[which(W_desc[1:k_hat] > 0)])
  } else{
    selected <- NULL
  }
  fdp_est <- ifelse(k_hat > 0, fdp[k_hat], 1)
  W_k_hat <- abs(W_desc[max(k_hat, 1)])

  results <- list(selected = selected, fdp_est = fdp_est,  W_k_hat = W_k_hat)
}

# determine if knockoff thinks j is promising in being non-null
# this functio may be used in screening (currently not).
kn_prefer <- function(j, kn_statistics, alpha, relax_factor = 1.5){
  p <- length(kn_statistics)
  prefer <- ((p+1-rank(abs(kn_statistics)))[j] <= 2/alpha * relax_factor)
  if(!prefer){
    kn_statistics[j] <- abs(kn_statistics[j])
    prefer <- (j %in% kn.select(kn_statistics, alpha * relax_factor,
                                selective = T, early_stop = T)$selected)
  }
  return(prefer)
}


