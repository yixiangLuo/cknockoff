# MRC_Xkn <- function(X, method = "mvr", normalize = T, num_processes = 1){
#   knockpy <- import("knockpy")
#
#   Simga_X <- t(X) %*% X
#   n <- NROW(X)
#   solver <- ifelse(method == "sdp", "sdp", "cd")
#
#   s_mat <- knockpy$smatrix$compute_smatrix(Simga_X, method = method, solver = solver,
#                                            num_processes = num_processes)
#   X_kn <- knockpy$knockoffs$produce_FX_knockoffs(X, Simga_X, s_mat) %>% matrix(nrow = n)
#
#   if(normalize){
#     X_kn <- scale(X_kn, center = FALSE, scale = sqrt(colSums(X_kn^2)))
#   }
#
#   return(X_kn)
# }

kn_rej_fdp <- function(kn_statistics, alpha, selective = TRUE,
                       k_min = 0, W_min = Inf, early_stop = F){
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
  k_hat <- max(k_min, ifelse(length(ok) > 0, max(ok), 0))

  if(W_min < Inf){
    ok <- which(abs(W_desc) >= W_min)
    k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
  }
  if(early_stop == 1){
    ok <- which(cumsum(W_desc > 0) < 1/alpha)
    k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
  }
  if(early_stop == 2){
    ok <- which((fdp <= 1) & (cumsum(W_desc > 0) < 1/alpha))
    k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
  }

  if(k_hat > 0){
    selected <- sort(W_desc_ind[which(W_desc[1:k_hat] > 0)])
  } else{
    selected <- NULL
  }
  fdp_est <- ifelse(k_hat > 0, fdp[k_hat], 1)
  W_k_hat <- abs(W_desc[max(k_hat, 1)])


  results <- list(selected = selected, fdp_est = fdp_est,  W_k_hat = W_k_hat)
}


