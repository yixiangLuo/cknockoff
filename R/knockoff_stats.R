#' Importance statistics based on Lasso penalty and Lagrange Multipliers
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y y vector of length n, containing the response variables.
#' @param sigma_tilde An estimator of the noise level independent of
#' \eqn{[X, X_k]^T y}.
#' By default, it is set to be NULL and will be computed inside the function
#' based on X, X_k, and y.
#'
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @details
#'
#' @family statistics
#'
#' @export
#'
#' @examples
stat.glmnet_lambdasmax_lm <- function(X, X_k, y, sigma_tilde = NULL) {
  n <- NROW(X)
  p <- NCOL(X)
  orig <- 1:p

  if(is.null(sigma_tilde)){
    sigma_tilde <- estimate_sigma(X, X_k, y)
  }
  lambda_min <- 2 * sigma_tilde / n

  # Randomly swap columns of X and Xk
  swap <- rbinom(ncol(X), 1, 0.5)
  swap_index <- c((1:p) + swap*p, ((p+1):(2*p)) - swap*p)
  X_full <- cbind(X, X_k)[, swap_index]

  # Compute statistics
  Z <- lasso_max_lambda_lm_glmnet(X_full, y, lambda_min = lambda_min)

  W <- pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W <- W * (1-2*swap)
}

lasso_max_lambda_lm_glmnet <- function(X, y, lambda_min, nlambda = 10) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)

  n <- NROW(X)
  p <- NCOL(X)

  lambda_max <- max(abs(t(X) %*% y)) / n
  if(lambda_min >= lambda_max){
    lambda_min <- lambda_max/1.2
  }
  nlambda <- 10
  k <- (0:(nlambda-1)) / nlambda
  lambda <- lambda_max * (lambda_min/lambda_max)^k


  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=T, standardize=F, standardize.response=F, family="gaussian")

  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL


  Z <- ifelse(is.na(indices), 0, fit$lambda[indices] * n)

  y_res <- y - (X %*% fit$beta + rep(fit$a0, each = n))
  # y_res <- y - predict(fit, X, lambda, type = "response")
  LM <- abs(t(X) %*% y_res)
  Z_cands <- c(lambda, 0) * n
  epsilon <- min(abs(diff(Z_cands)), 1) * 1e-6
  for(Z_i in (2:length(Z_cands))){
    this_level <- which(abs(Z - Z_cands[Z_i]) < epsilon)
    LM_before <- LM[this_level, Z_i - 1]

    exceed_ind <- which(LM_before > Z_cands[Z_i-1] - epsilon)
    exceed <- LM_before[exceed_ind] - (Z_cands[Z_i-1] - epsilon)
    if(length(exceed_ind) > 0){
      LM_before[exceed_ind] <- (Z_cands[Z_i-1] - epsilon) + exceed * epsilon / max(abs(exceed))
    }
    exceed_ind <- which(LM_before < Z_cands[Z_i] + epsilon)
    exceed <- LM_before[exceed_ind] - (Z_cands[Z_i] + epsilon)
    if(length(exceed_ind) > 0){
      LM_before[exceed_ind] <- (Z_cands[Z_i] + epsilon) + exceed * epsilon / max(abs(exceed))
    }

    Z[this_level] <- LM_before
  }

  return(Z)
}


#' Importance statistics based on Lasso coefficients and Lagrange Multipliers
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y y vector of length n, containing the response variables.
#' @param sigma_tilde An estimator of the noise level independent of
#' \eqn{[X, X_k]^T y}.
#' By default, it is set to be NULL and will be computed inside the function
#' based on X, X_k, and y.
#'
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @details
#'
#' @family statistics
#'
#' @export
#'
#' @examples
stat.glmnet_coefdiff_lm <- function(X, X_k, y, sigma_tilde = NULL) {
  n <- NROW(X)
  p <- NCOL(X)
  orig <- 1:p

  if(is.null(sigma_tilde)){
    sigma_tilde <- estimate_sigma(X, X_k, y)
  }
  lambda <- 2 * sigma_tilde / n

  # Randomly swap columns of X and Xk
  swap <- rbinom(ncol(X), 1, 0.5)
  swap_index <- c((1:p) + swap*p, ((p+1):(2*p)) - swap*p)
  X_full <- cbind(X, X_k)[, swap_index]

  # Compute statistics
  fit <- glmnet::glmnet(X_full, y, lambda = lambda,
                        intercept = T, standardize = F,
                        standardize.response = F, family = "gaussian")

  Z <- abs(c(as.matrix(fit$beta)))
  W <- Z[orig] - Z[orig+p]

  y_res <- as.matrix(y - (X_full %*% fit$beta + rep(fit$a0, each = n)))
  LM <- c(abs(t(y_res) %*% X_full))
  LM_add_on <- rep(NA, p)
  not_in_model <- which(W == 0)
  LM_add_on[not_in_model] <- LM[not_in_model] - LM[not_in_model+p]
  in_model <- which(W != 0)
  LM_add_on[in_model] <- lambda * n * sign(W[in_model])

  W <- W + LM_add_on

  # Correct for swapping of columns of X and Xk
  W <- W * (1-2*swap)
}

estimate_sigma <- function(X, X_k, y){
  n <- NROW(X)
  p <- NCOL(X)

  if(n >= 2*p + 1){
    QR <- qr(cbind(X, X_k))
    Q <- qr.Q(QR, complete = F)
    sigma_tilde <- sqrt((sum(y^2) - sum((matrix(y, nrow = 1) %*% Q)^2)) / (n - 2*p))
  } else if(n == 2*p){
    QR <- qr(X)
    Q <- qr.Q(QR, complete = F)
    sigma_tilde <- sqrt((sum(y^2) - sum((matrix(y, nrow = 1) %*% Q)^2)) / (n - p))
  } else{
    stop("X must have dimensions n >= 2*p.")
  }

  return(sigma_tilde)
}
