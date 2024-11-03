#' Lasso Coefficient-Difference feature statistics with tiebreaker
#'
#' Fit Lasso on the augmented model \eqn{y ~ [X, \widetilde{X}]} with
#' regularization parameter \eqn{\lambda}.
#' Then, compute the difference statistic
#'   \deqn{W^0_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
#' jth variable and its knockoff, respectively.
#' For those variables that both themselves and their knockoffs are not selected
#' by Lasso, we break the tie by their correlation with the residue, by defining
#'   \deqn{\rho_j = |X_j^T residue| - |\tilde{X}_j^T residue|}
#' for them and \eqn{\rho_j = 2 \lambda} for the others.
#' The final feature statistics are
#' \deqn{
#'   W_j = W^0_j + \rho_j.
#' }
#' \eqn{\lambda} is set to be \eqn{2 \tilde{\sigma}}, where \eqn{\tilde{\sigma}}
#' is an estimator of the noise level that is independent of
#' \eqn{[X, \widetilde{X}]^T y}.
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y y vector of length n, containing the response variables.
#' @param sigma_tilde An estimator of the noise level that is independent of
#' \eqn{[X, \widetilde{X}]^T y}.
#' If not provided, it will be computed inside the function
#' based on X, X_k, and y, in which case we must have \eqn{n > 2p}.
#'
#' @return A vector of knockoff feature statistics \eqn{W} of length p.
#'
#' @details
#'
#' If sigma_tilde is not provided and \eqn{n = 2p}, sigma_tilde will be computed
#' from the residue of the OLS fitting \eqn{y ~ X}.
#' The resulting estimator is thus not independent of \eqn{[X, \widetilde{X}]^T y}.
#' Users should avoid this case if they want a guaranteed FDR control.
#' Though in practice it shouldn't really make a difference.
#'
#' Details of the calculation of this feature statistics can be found in (our paper).
#'
#' The implementation of this function is modified from the
#' \code{knockoff::stat.glmnet_coefdiff()} function.
#'
#' @family statistics
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 2.5 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#' print(which(1:p %in% nonzero))
#'
#' result <- cknockoff(X, y,
#'                     knockoffs = ckn.create.fixed,
#'                     statistic = stat.glmnet_coefdiff_tiebreak,
#'                     alpha = 0.05, n_cores = 1)
#' print(result$selected)
#'
#' @export
stat.glmnet_coefdiff_tiebreak <- function(X, X_k, y, sigma_tilde = NULL) {
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
                        intercept = F, standardize = F,
                        standardize.response = F, family = "gaussian")

  Z <- abs(c(as.matrix(fit$beta)))
  W <- Z[orig] - Z[orig+p]

  y_res <- y - as.matrix(X_full %*% fit$beta) + rep(fit$a0, each = n)
  res_corr <- c(abs(t(y_res) %*% X_full))
  tie_breaker <- lambda * n * sign(W)
  not_in_model <- which(W == 0)
  tie_breaker[not_in_model] <- res_corr[not_in_model] - res_corr[not_in_model+p]

  W <- W + tie_breaker

  # Correct for swapping of columns of X and Xk
  W <- W * (1-2*swap)
}


#' Efficient Lasso Signed-Max feature statistics
#'
#' This function provides a computationally efficient feature statistic that
#' behave similarly to the Lasso Signed-Max feature statistics produced by
#' the \code{knockoff::stat.glmnet_lambdasmax()} function.
#' We compute a \eqn{\hat{\lambda}_j} as an approximation of the maximal
#' regularization parameter \eqn{\lambda} at which the jth variable enters the
#' Lasso model.
#' And obtain the feature statistics
#' \deqn{
#' W_j = (\hat{\lambda}_{j} \vee \hat{\lambda}_{j+m}) \cdot
#' sgn (\hat{\lambda}_{j} - \hat{\lambda}_{j+m})
#' }
#'
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y y vector of length n, containing the response variables.
#' @param sigma_tilde An estimator of the noise level that is independent of
#' \eqn{[X, \widetilde{X}]^T y}.
#' If not provided, it will be computed inside the function
#' based on X, X_k, and y, in which case we must have \eqn{n > 2p}.
#' @param nlambda the number of grid points in computing the Lasso path.
#' A larger value will make the calculation more precise but more expensive.
#' The default value is 10.
#'
#' @return A vector of knockoff feature statistics \eqn{W} of length p.
#'
#' @details
#'
#' If sigma_tilde is not provided and \eqn{n = 2p}, sigma_tilde will be computed
#' from the residue of the OLS fitting \eqn{y ~ X}.
#' The resulting estimator is thus not independent of \eqn{[X, \widetilde{X}]^T y}.
#' Users should avoid this case if they want a guaranteed FDR control.
#' Though in practice it shouldn't really make a difference.
#'
#' Details of the calculation of this feature statistics can be found in (our paper).
#'
#' The implementation of this function is modified from the
#' \code{knockoff::stat.glmnet_lambdasmax()} function.
#'
#' @family statistics
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 2.5 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#' print(which(1:p %in% nonzero))
#'
#' result <- cknockoff(X, y,
#'                     knockoffs = ckn.create.fixed,
#'                     statistic = stat.glmnet_lambdasmax_coarse,
#'                     alpha = 0.05, n_cores = 1)
#' print(result$selected)
#'
#' @export
stat.glmnet_lambdasmax_coarse <- function(X, X_k, y,
                                          sigma_tilde = NULL,
                                          nlambda = 10) {
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
  Z <- lasso_max_lambda_lm_glmnet(X_full, y,
                                  lambda_min = lambda_min,
                                  nlambda = nlambda)

  W <- pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W <- W * (1-2*swap)
}

lasso_max_lambda_lm_glmnet <- function(X, y, lambda_min, nlambda = 10) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)

  n <- NROW(X)
  p <- NCOL(X)

  lambda_max <- max(abs(t(X) %*% y))
  lambda_min <- min(lambda_min, lambda_max/2)

  eta <- (lambda_min/lambda_max)^(1/(nlambda-1))
  lambda_grid <- lambda_max * eta^(0:(nlambda-1))

  fit <- glmnet::glmnet(X, y, lambda = lambda_grid / n,
                        intercept = F, standardize = F, standardize.response = F,
                        family="gaussian")

  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  time_of_entry <- apply(fit$beta, 1, first_nonzero)
  names(time_of_entry) <- NULL
  time_of_entry[is.na(time_of_entry)] <- nlambda + 1

  y_res <- c(y) - as.matrix(X %*% fit$beta) + rep(fit$a0, each = n)
  tie_breaker <- sapply(1:p, function(j){
    if(time_of_entry[j] > 1){
      return(abs(t(X[, j]) %*% y_res[, time_of_entry[j]-1]))
    } else{
      return(0)
    }
  })

  iota <- min(abs(diff(lambda_grid)), 1) * 1e-6
  lambda_plus <- c(lambda_max / eta, lambda_grid)[time_of_entry]
  lambda_enter <- c(lambda_grid, 0)[time_of_entry]
  lambda_hat <- pmax(tie_breaker - iota, lambda_enter) + iota * tie_breaker / lambda_plus
  # lambda_hat <- lambda_plus * eta^(1 - tie_breaker / lambda_plus)

  return(lambda_hat)
}




# estimate sigma from the residue of regressing y ~ [X, X_k]
estimate_sigma <- function(X, X_k, y){
  n <- NROW(X)
  p <- NCOL(X)

  intercept <- (max(abs(c(colMeans(X), colMeans(X_k)))) < 1e-6)
  if(intercept) y <- scale(y, center = T, scale = F)

  if(n >= 2*p + 1 + intercept){
    QR <- qr(cbind(X, X_k))
    Q <- qr.Q(QR, complete = F)
    sigma_tilde <- sqrt((sum(y^2) - sum((matrix(y, nrow = 1) %*% Q)^2)) / (n - 2*p - intercept))
  } else if(n <= 2*p + intercept){
    QR <- qr(X)
    Q <- qr.Q(QR, complete = F)
    sigma_tilde <- sqrt((sum(y^2) - sum((matrix(y, nrow = 1) %*% Q)^2)) / (n - p - intercept))
  } else{
    if(intercept){
      stop("X must have dimensions n >= 2*p+1 when colMeans([X, Xk])=0 (intercept included effectively).")
    } else{
      stop("X must have dimensions n >= 2*p.")
    }
  }

  return(sigma_tilde)
}
