parse_args <- function(X, y, knockoffs, statistic,
                       alpha, n_cores, knockoff.type,
                       prelim_result, X.pack, envir){

  if(is(prelim_result, "cknockoff.result")){
    X <- prelim_result$X
    y <- prelim_result$y
    knockoffs <- prelim_result$Xk
    statistic <- prelim_result$record$statistic
    alpha <- prelim_result$record$alpha
    X.pack <- prelim_result$record$X.pack
    record <- c(prelim_result$record[c("iteration", "kn.selected",
                                       "checked_so_far", "next_check_num")],
                list(kn.statistic = prelim_result$kn.statistic,
                     selected_so_far = prelim_result$selected))

    if(evalq(missing(n_cores), parent.frame())){
      n_cores <- prelim_result$record$n_cores
    }
  }
  else if(is(prelim_result, "knockoff.result")){
    X <- prelim_result$X
    y <- prelim_result$y
    knockoffs <- prelim_result$Xk
    # X.pack <- NULL
    record <- list(iteration = 0,
                   kn.statistic = prelim_result$statistic,
                   kn.selected = prelim_result$selected)

    statistic_missing <- evalq(missing(statistic), parent.frame())
    tryCatch(
      {
        statistic <- eval(prelim_result$call$statistic, envir = envir)
      },
      error = function(error_message){
        if(statistic_missing){
          stop("Cannot read function \"statistic\" from the knockoff.result object \"prelim_result\". \n",
               "please supply it explicitly by passing \"statistic = ...\" in cknockoff().\n",
               "Please note you must use the same \"statistic\" function as in knockoff to get valid inference.")
        }
        # else use arg statistic
      }
    )

    alpha_missing <- evalq(missing(alpha), parent.frame())
    tryCatch(
      {
        alpha <- eval(prelim_result$call$fdr, envir = envir)
      },
      error = function(error_message){
        if(alpha_missing){
          stop("Cannot read nominated \"fdr\" from the knockoff.result object \"prelim_result\". \n",
               "please supply it explicitly by passing \"alpha = ...\" in cknockoff().\n",
               "Please note you must use the same nominated FDR level as in knockoff to get valid inference.")
        }
        # else use arg alpha
      }
    )
  }
  else if((!missing(X) || !is.null(X.pack)) && !missing(y)){
    record <- NULL
  }
  else{
    stop("X/X.pack and y must be provided if knockoff/cknockoff.result object isn't.")
  }

  n_cores <- as.integer(n_cores)

  # for validation
  if(is(X.pack, "cknockoff.X.pack")){
    X <- X.pack$X
    knockoffs <- X.pack$X_kn
  }

  return(list(X = X, y = y, knockoffs = knockoffs, statistic = statistic,
              alpha = alpha, n_cores = n_cores, X.pack = X.pack,
              record = record))
}

check_args <- function(args){
  # Validate input types
  if(!is.data.frame(args$X) && !is.matrix(args$X)){
    stop("X must be a numeric matrix or data frame.")
  }
  if(!is.numeric(as.matrix(args$X))){
    stop("X must be a numeric matrix or data frame.")
  }

  if(!is.numeric(args$y)){
    stop("y must be of numeric type.")
  }

  if(!is.function(args$knockoffs) && !is.matrix(args$knockoffs)){
    stop("Input knockoffs must be either a function or matrix.")
  }
  if(!is.function(args$statistic)){
    stop("Input statistic must be a function.")
  }

  if(!is.numeric(args$alpha) || (args$alpha < 0 || args$alpha > 1)){
    stop("alpha must be a number between 0 and 1.")
  }

  if(!is.integer(args$n_cores) || args$n_cores <= 0){
    stop("Input n_cores must be a positive integer.")
  }

  if(!is.null(args$X.pack) && !is(args$X.pack, "cknockoff.X.pack")){
    stop("X.pack must be a cknockoff.X.pack object if provided.")
  }

  # Validate input dimensions

  y.length <- length(c(args$y))
  if(is(args$X.pack, "cknockoff.X.pack")){
    if(y.length != NROW(args$X.pack$X) && y.length != args$X.pack$X.org.nrow){
      stop("X and y have inconsistent number of rows.")
    }
  } else if(y.length != NROW(args$X)){
    stop("X and y have inconsistent number of rows.")
  }

  invisible()
}

process_args <- function(args){
  if(!is(args$X.pack, "cknockoff.X.pack")){
    args$X.pack <- process_X(args$X, knockoffs = args$knockoffs)
  }

  args$y.pack <- process_y(args$X.pack, args$y, randomize = F)

  # # not a good practice. Now this is a infinite self-reference.
  # # even after fixing this, arg is a large variable live in this local environment
  # # that can't be removed.
  # args_name <- names(formals(args$statistic))
  # if(!("sigma_tilde" %in% args_name)){
  #   args$statistic <- function(X, X_k, y, sigma_tilde){
  #     args$statistic(X, X_k, y)
  #   }
  # }

  args$parallel <- process_parallel(args$n_cores)

  return(args)
}


#' Compute and store useful matrices based on X
#'
#' @param X X n-by-p matrix of original variables.
#' @param knockoffs either knockoff matrix of X or a knockoffs function that can
#' generate it.
#' If the knockoff matrix is supplied, both X and knockoffs should be properly
#' normalized, e.g. using the returned X and Xk of the knockoff::create.fixed function.
#' By default, knockoff::create.fixed is used.
#'
#' @return An object of class "cknockoff.X.pack". This object is a list containing
#' many matrices like the knockoff matrix and a basis of the linear space spanned by X, etc.
#' The users don't need to work on the object themselves but pass it to the "X.pack"
#' parameter in cknockoff() if needed.
#' @export
#'
#' @examples
process_X <- function(X, knockoffs = knockoff::create.fixed){
  if(NCOL(X) <= 1){
    stop("X must have at least two columns.")
  }
  if(NROW(X) <= NCOL(X)){
    stop("X must have dimensions n > p")
  }

  if(is.data.frame(X)){
    X.names <- names(X)
    X <- as.matrix(X, rownames.force = F)
  } else if(is.matrix(X)){
    X.names <- colnames(X)
  } else{
    stop("X must be a numeric matrix or data frame.")
  }

  X.org.nrow <- NROW(X)
  p <- NCOL(X)

  if(is.matrix(knockoffs)){
    # when X is raw but X_kn is augmented
    if(!identical(dim(X), dim(knockoffs))){
      stop("X and its knockoff matrix (knockoffs) have inconsistent dimensions. ",
           "Is the knockoff matrix augmented?")
    }
    # when X is raw but X_kn is scaled by create.fixed
    if(abs(sum(X[, 1] * X[, 2]) - sum(knockoffs[, 1] * knockoffs[, 2])) > 1e-5){
      stop("the provided knockoffs is not a valid knockoff matrix of X. ",
           "Do they have consistent normalization?")
    }
    X_kn <- knockoffs
  } else if(is.function(knockoffs)){
    # augment X beforehead, otherwise y or sigma is needed in create.fixed
    if(X.org.nrow < 2*p){
      X <- rbind(X, matrix(0, 2*p-X.org.nrow, p))
    }
    kn_variables <- knockoffs(X)
    X <- kn_variables$X
    X_kn <- kn_variables$Xk
  } else{
    stop("knockoff matrix/generating method of incorrect type.")
  }

  n <- NROW(X)

  # augment 2*p+1 in addition to 2*p to make sure X_kn has the last row = 0
  if(n < 2*p+1){
    warning('Input X has dimensions n < 2p+1. ',
            'Augmenting the model with extra rows.', immediate.=T)
    X <- rbind(X, matrix(0, 2*p+1-n, p))
    X_kn <- rbind(X_kn, matrix(0, 2*p+1-n, p))
    n <- 2*p+1
  }

  # if(method == "sdp"){
  #   X_kn <- knockoff::create.fixed(X, method = "sdp")$Xk
  # } else{
  #   X_kn <- MRC_Xkn(X, method = method, normalize = T)
  # }

  QR <- qr(cbind(X, X_kn))

  pivot_back <- 1:p
  pivot_back[QR$pivot] <- pivot_back

  Q_XXk <- qr.Q(QR, complete = T)
  R_XXk <- qr.R(QR)

  if(min(abs(diag(R_XXk)[pivot_back[1:p]])) < 1e-7){
    stop("X doesn't have full column rank.")
  }

  if(!identical(QR$pivot[1:p], 1:p)){
    X.pack.data <- process_X_robust(X, X_kn)
  } else{
    X_basis <- Q_XXk[, 1:p]
    X_res_Xk_basis <- Q_XXk[, (p+1):(2*p)]
    XXk_res_unit <- Q_XXk[, (2*p+1)]

    # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
    vj_mat <- sapply(1:p, function(j){
      # Q_X[, j] = X_j orthogonal projected onto X_{1:j-1}
      # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
      #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
      # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
      coord <- forwardsolve(t(R_XXk[j:p,j:p]), c(1, rep(0, p-j)))
      vj <- Q_XXk[, j:p] %*% matrix(coord, nrow = p-j+1)
      vj <- vj / sqrt(sum(vj^2))
    })

    X.pack.data <- list(X = X, X_kn = X_kn,
                        X_basis = X_basis, vj_mat = vj_mat,
                        X_res_Xk_basis = X_res_Xk_basis,
                        XXk_res_unit = XXk_res_unit)

  }

  X.pack <- structure(c(X.pack.data, list(X.names = X.names, X.org.nrow = X.org.nrow)),
                       class = "cknockoff.X.pack")

  return(X.pack)

}

process_X_robust <- function(X, X_kn){
  n <- NROW(X)
  p <- NCOL(X)

  QR_X <- qr(X)

  Q_X <- qr.Q(QR_X, complete = T)
  R_X <- qr.R(QR_X)

  # qr(X) is possibly pivoted. The following matrices are only used as bases for the subspace.
  X_basis <- Q_X[, 1:p]

  # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
  vj_mat <- sapply(1:p, function(j){
    # Q_X[, j] = X_j orthogonal projected onto X_{1:j-1}
    # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
    #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
    # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
    coord <- forwardsolve(t(R_X[j:p,j:p]), c(1, rep(0, p-j)))
    vj <- Q_X[, j:p] %*% matrix(coord, nrow = p-j+1)
    vj <- vj / sqrt(sum(vj^2))
  })
  # pivot back
  vj_mat[, QR_X$pivot] <- vj_mat

  X_res_Xk <- X_kn - X_basis %*% (t(X_basis) %*% X_kn)

  # run QR for X_res_Xk to avoid pivoting between X and Xk, which makes the calculation of vj_mat invalid.
  QR_Xk <- qr(X_res_Xk)
  Q_Xk <- qr.Q(QR_Xk, complete = F)
  X_res_Xk_basis <- Q_Xk[, 1:p]

  Q_XXk <- qr.Q(qr(cbind(X_basis, X_res_Xk_basis)), complete = T)
  XXk_res_unit <- Q_XXk[, 2*p+1]

  X.pack.data <- list(X = X, X_kn = X_kn,
                      X_basis = X_basis, vj_mat = vj_mat,
                      X_res_Xk_basis = X_res_Xk_basis,
                      XXk_res_unit = XXk_res_unit)

  return(X.pack.data)
}

process_y <- function(X.pack, y, randomize = F){

  n <- NROW(X.pack$X)
  p <- NCOL(X.pack$X)
  df <- n-p

  y.org.nrow <- length(y)

  y <- matrix(y, nrow = 1)
  y_Pi_X <- X.pack$X_basis %*% t(y %*% X.pack$X_basis[1:y.org.nrow, ])
  y_Pi_X_res_norm2 <- sum(y^2) - sum(y_Pi_X^2)

  if(y.org.nrow < 2*p){
    if(randomize){
      y.extra <- rnorm(2*p-y.org.nrow,
                       sd = sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p)))
    } else{
      y.extra <- with_seed(0, rnorm(2*p-y.org.nrow,
                                    sd = sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p))))
    }
    y <- c(y, y.extra)
  }
  if(y.org.nrow < 2*p+1){
    y <- c(y, sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p)))
  }
  y_Pi_X_res_norm2 <- sum(y^2) - sum(y_Pi_X^2)

  y <- matrix(y, nrow = 1)


  vjy_obs <- c(y %*% X.pack$vj_mat)

  y_Pi_Xnoj <- sapply(1:p, function(j){
    y_Pi_X - X.pack$vj_mat[, j] * vjy_obs[j]
  })

  y_Pi_Xnoj_res_norm2 <- y_Pi_X_res_norm2 + vjy_obs^2
  sigmahat_XXk_res <- sqrt((y_Pi_X_res_norm2 - sum((y %*% X.pack$X_res_Xk_basis)^2)) / (n - 2*p))

  Xy_bias <- rep(NA, p)

  return(list(n = n, p = p, df = df, y = as.vector(y), vjy_obs = vjy_obs,
              y_Pi_X_res_norm2 = y_Pi_X_res_norm2, y_Pi_Xnoj_res_norm2 = y_Pi_Xnoj_res_norm2,
              sigmahat_XXk_res = sigmahat_XXk_res,
              y_Pi_Xnoj = y_Pi_Xnoj,
              Xy_bias = Xy_bias))
}

process_parallel <- function(n_cores){
  if(n_cores > 1 && !requireNamespace("doParallel", quietly=T)) {
    warning("doParallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("parallel", quietly=T)) {
    warning("parallel is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1 && !requireNamespace("foreach", quietly=T)) {
    warning("foreach is not installed. Using sequential computing instead.", call. = F, immediate. = T)
    n_cores <- 1
  }
  if(n_cores > 1){
    cores_avail <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if(n_cores > cores_avail){
      warning(paste("The requested number of cores is not available. Using instead", cores_avail, "cores."), immediate. = T)
      n_cores <- cores_avail
    }
  }
  if(n_cores > 1){
    doParallel::registerDoParallel(n_cores)

    forall <- foreach::foreach
    `%exec%` <- foreach::`%dopar%`
  } else{
    forall <- iterate_seq
    `%exec%` <- `%do_seq%`
  }

  return(list(iterator = forall, connector = `%exec%`))
}

