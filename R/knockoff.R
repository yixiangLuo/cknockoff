#' Title
#'
#' @param X
#' @param method
#' @param solver
#' @param knockpy
#' @param num_processes
#'
#' @return a list of (X, X_kn)
#' @export
create.fixed.MRC <- function(X,
                             method = c("mvr", "maxent", "mmi",
                                        "sdp", "equicorrelated", "ci"),
                             solver = NULL,
                             knockpy = NULL,
                             num_processes = 1){

  if(!requireNamespace("reticulate", quietly=T)) {
    stop("reticulate is required for calling python package 'knockpy' but is not installed.")
  }

  method <- match.arg(method)
  if(is.null(solver)){
    solver <- ifelse(method == "sdp", "sdp", "cd")
  }

  if(is.null(knockpy)){
    knockpy <- reticulate::import("knockpy")
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

#' Title
#'
#' @param kn_statistics Knockoff W-statistics
#' @param alpha nominated FDR level
#' @param selective T: use selective SeqStep, F: use SeqStep
#' @param early_stop T/F
#'
#' @return list(selected, estimated fdp, |W_(k hat)|)
#' @export
#'
#' @examples
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


