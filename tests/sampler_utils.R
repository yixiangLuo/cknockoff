# naively sampling y conditional on Sj
y_condj_sample <- function(y, X, j, sample_size, seed = 1, q_cap = Inf){
  set.seed(seed + 2000)

  n <- NROW(X)

  projj_y <- lm(formula = y ~ X[, -j] + 0)$fitted.values
  radius <- sqrt(sum(y^2) - sum(projj_y^2))

  y_sample <- matrix(rnorm(n * sample_size), nrow = n)

  y_sample <- y_sample - lm(formula = y_sample ~ X[, -j] + 0)$fitted.values
  y_sample <- scale(y_sample, center = FALSE, scale = sqrt(colSums(y_sample^2)) / radius)
  y_sample <- projj_y + y_sample

  return(y_sample)
}


# generate X matrix
gene_X <- function(X_type = "IID_Normal", n, p, X_seed = 1){
  cor_radius <- 5
  set.seed(X_seed)
  if(X_type == "IID_Normal"){
    X <- matrix(rnorm(n*p), n) / sqrt(n)
  } else if(X_type == "Coef_AR"){
    rho <- 0.5

    cov_mat <- solve(rho^(abs(outer(1:p, 1:p, "-"))))

    R <- chol(cov_mat)
    basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
    X <- basis %*% R
  } else if(X_type == "X_AR"){
    rho <- 0.5

    cov_mat <- rho^(abs(outer(1:p, 1:p, "-")))

    R <- chol(cov_mat)
    basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
    X <- basis %*% R
  } else if(X_type == "Homo_Block"){
    rho <- 0.5
    block_size <- 10

    blockSigma <- matrix(rho, block_size, block_size)
    diag(blockSigma) <- 1

    cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)

    R <- chol(cov_mat)
    basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
    X <- basis %*% R
  } else if(X_type == "MCC"){
    if(n %% (p+1) == 0){
      X <- lapply(1:(n/(p+1)), function(i){
        rbind(diag(rep(1, p)), rep(0, p))
      })
      X <- do.call(rbind, X)
      X <- scale(X, center = T, scale = F)
      X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
    } else{
      cov_mat <- matrix(-1/p, p, p)
      diag(cov_mat) <- 1

      R <- chol(cov_mat)
      basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
      X <- basis %*% R
    }
  } else if(X_type == "MCC_Block"){
    block_size <- 5

    blockSigma <- matrix(-1/block_size, block_size, block_size)
    diag(blockSigma) <- 1

    cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)

    R <- chol(cov_mat)
    basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
    X <- basis %*% R
  } else if(X_type == "Sparse"){
    X <- diag(1, nrow = n, ncol = p)
    nonzeros <- matrix(NA, nrow = p, ncol = 2)
    nonzeros[, 1] <- 1:p
    nonzeros[, 2] <- sample(1:p, p, replace = T)
    X[nonzeros] <- X[nonzeros] + 0.5
  }
  # X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))

  return(X)
}

# generate correlated noise
corr_noise <- function(n, rho){
  cov_mat <- matrix(rho, n, n)
  diag(cov_mat) <- 1
  # cov_mat <- matrix(0, n, n)
  # cov_mat[abs(row(cov_mat) - col(cov_mat)) <= 1] <- -rho * 0.5
  # diag(cov_mat) <- 1

  R <- chol(cov_mat)
  noise <- t(R) %*% matrix(rnorm(n), nrow = n)

  return(noise)
}


# compute the t-vals of a linear regression test problem
lm_to_t <- function(y, X, Sigma = NULL){
  n <- NROW(X)
  p <- NCOL(X)

  if(is.null(Sigma)){
    Sigma <- solve(t(X) %*% X)
  }
  Xy <- t(X) %*% y
  df <- n - p

  zvals <- Sigma %*% Xy
  sigmahat <- as.numeric(sqrt((sum(y^2) - t(Xy) %*% zvals) / df))
  tvals <- zvals / sqrt(diag(Sigma)) / sigmahat

  return(list(tvals = tvals, df = df))
}

# apply BH method to a linear regression test problem
BH_lm <- function(y, X, side = "two", alpha,
                  weights = rep(1, NCOL(X)) / NCOL(X),
                  Sigma = NULL){

  t_result <- lm_to_t(y, X, Sigma)
  pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")

  BH_result <-  BH_weighted(pvals, alpha, weights)

  return(BH_result)
}

## calibrate signal strength, modified from Lihua
BH_lm_calib <- function(X, pi1, noise = quote(rnorm(n)),
                        mu_posit_type, mu_size_type,
                        side,
                        nreps = 1000,
                        alpha = 0.05,
                        target = 0.3,
                        n_cores = 7){
  n <- nrow(X)
  p <- ncol(X)
  Sigma <- solve(t(X) %*% X)
  H <- X %*% Sigma %*% t(X)
  df <- n - p

  beta_list <- lapply(1:nreps, function(i){
    beta <- genmu(p, pi1, 1, mu_posit_type, mu_size_type)
    if (side == "right"){
      beta <- abs(beta)
    } else if (side == "left"){
      beta <- -abs(beta)
    }
    return(beta)
  })
  eps_list <- lapply(1:nreps, function(i){
    eval(noise)
  })

  registerDoParallel(n_cores)

  BH_power <- function(mu1){
    power <- unlist(foreach(i = 1:nreps) %dopar% {
      H0 <- beta_list[[i]] == 0
      beta <- beta_list[[i]] * mu1
      eps <- eps_list[[i]]
      y <- X %*% beta + eps

      rejs_BH <- BH_lm(y, X, side = "two", alpha, Sigma = Sigma)$rejs
      power_sample <- calc_FDP_power(rejs_BH, H0)[2]

      return(power_sample)
    })
    mean(power) - target
  }

  lower <- 0
  upper <- 10
  while (TRUE & upper < 1000){
    tmp <- try(uniroot(BH_power, c(lower, upper))$root)
    if (class(tmp) == "try-error"){
      upper <- upper * 2
    } else {
      return(tmp)
    }
  }
  return(NA)
}

# generate beta, borrowed from Lihua
genmu <- function(n, pi1, mu1,
                  posit_type = c("random", "fix"),
                  mu_type = 1:3){
  m <- ceiling(n * pi1)
  posit_type <- posit_type[1]
  mu_type <- mu_type[1]
  if (posit_type == "random"){
    inds <- seq(1, n, floor(1 / pi1))[1:m]
  } else if (posit_type == "fix"){
    inds <- 1:m
  }
  mu <- rep(0, n)
  altmu <- switch(mu_type,
                  `1` = rep(1, m),
                  `2` = rnorm(m),
                  `3` = rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1))
  mu[inds] <- mu1 * altmu
  mu
}

# compute FDP and power
calc_FDP_power <- function(rejs, H0, sign_predict = NULL, sign_beta = NULL){
  nrejs <- length(rejs)

  false_discovery <- length(intersect(rejs, which(H0)))
  true_discovery <- nrejs - false_discovery

  FDP <- false_discovery / max(nrejs, 1)
  power <- true_discovery / max(sum(!H0), 1)

  if(!is.null(sign_predict)){
    false_dir <- sum((sign_predict * sign_beta)[rejs] <= 0)
    true_dir <- nrejs - false_dir

    FDP_dir <- false_dir / max(nrejs, 1)
    power_dir <- true_dir / length(sign_beta)
  } else{
    FDP_dir <- NA
    power_dir <- NA
  }

  return(c(FDP, power, FDP_dir, power_dir))
}
