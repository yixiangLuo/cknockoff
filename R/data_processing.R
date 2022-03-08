# parse the arguments supplied to cknockoff()
parse_args <- function(X, y, intercept,
                       knockoffs, statistic,
                       alpha, Rstar_refine,
                       n_cores,
                       prelim_result, X.pack,
                       envir){
  # when cknockoff.result is supplied
  if(is(prelim_result, "cknockoff.result")){
    X.pack <- prelim_result$record$X.pack
    y.pack <- prelim_result$record$y.pack
    intercept <- prelim_result$intercept
    X <- prelim_result$X
    y <- prelim_result$y
    knockoffs <- prelim_result$Xk
    statistic <- prelim_result$record$statistic
    alpha <- prelim_result$record$alpha

    record <- c(prelim_result$record[c("iteration", "kn.selected",
                                       "checked_so_far", "next_check_num")],
                list(kn.statistic = prelim_result$kn.statistic,
                     selected_so_far = prelim_result$selected))

    # make Rstar_refine=T if previously set TRUE, else make the non-rejected unchecked
    if(prelim_result$record$Rstar_refine == T) Rstar_refine <- T
    else if(Rstar_refine == T) record$checked_so_far <- prelim_result$selected

    if(evalq(missing(n_cores), parent.frame())){
      n_cores <- prelim_result$record$n_cores
    }
  }
  # when knockoff.result is supplied
  else if(is(prelim_result, "knockoff.result")){
    X <- prelim_result$X
    y <- prelim_result$y
    knockoffs <- prelim_result$Xk
    intercept <- (max(abs(c(colMeans(X), colMeans(knockoffs)))) < 1e-6)
    # X.pack <- NULL
    y.pack <- NA
    record <- list(iteration = 0,
                   kn.statistic = prelim_result$statistic,
                   kn.selected = prelim_result$selected)

    # try fetch the values of statistic and alpha by their names
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
    y.pack <- NA
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

  return(list(X = X, y = y, intercept = intercept,
              X.pack = X.pack, y.pack = y.pack,
              knockoffs = knockoffs, statistic = statistic,
              alpha = alpha, n_cores = n_cores,
              Rstar_refine = Rstar_refine, record = record))
}

# check if the arguments are valid
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
  if(is.na(args$y.pack)){
    y.length <- length(c(args$y))
    if(is(args$X.pack, "cknockoff.X.pack")){
      if(y.length != NROW(args$X.pack$X.org) &&
         y.length != args$X.pack$X.org.nrow){
        stop("X and y have inconsistent number of rows.")
      }
    } else if(y.length != NROW(args$X)){
      stop("X and y have inconsistent number of rows.")
    }
  }

  invisible()
}

# preprocess the X and y for better efficiency
process_args <- function(args){
  if(!is(args$X.pack, "cknockoff.X.pack")){
    args$X.pack <- process_X(args$X, knockoffs = args$knockoffs,
                             intercept = args$intercept)
  }
  if(is.na(args$y.pack)){
    y.data <- transform_y(args$X.pack, args$y,
                          intercept = args$intercept,
                          randomize = F)
    args$y.pack <- process_y(args$X.pack, y.data)
  }

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


#' Compute and store matrices needed by cknockoff() from on X
#'
#' @param X X n-by-p matrix of original variables.
#' @param knockoffs either knockoff matrix of X or a knockoffs function that can
#' generate it.
#' If the knockoff matrix is supplied, both X and knockoffs should be properly
#' normalized, e.g. using the returned X and Xk of the ckn.create.fixed function.
#' By default, ckn.create.fixed is used.
#'
#' @return An object of class "cknockoff.X.pack". This object is a list containing
#' many matrices like the knockoff matrix and a basis of the linear space spanned by X, etc.
#' Users don't need to work on the object themselves but pass it to the "X.pack"
#' parameter in cknockoff() if needed.
#'
#' @examples
#' p <- 100; n <- 300; k <- 15
#' X <- matrix(rnorm(n*p), n)
#' nonzero <- sample(p, k)
#' beta <- 2.5 * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#' print(which(1:p %in% nonzero))
#'
#' X.pack <- process_X(X, knockoffs = ckn.create.fixed)
#'
#' result <- cknockoff(X, y,
#'                     alpha = 0.05,
#'                     n_cores = 1,
#'                     X.pack = X.pack)
#' print(result$selected)
#'
#'
#' @export
process_X <- function(X, knockoffs = ckn.create.fixed,
                      intercept){
  # validate the arguments
  if(NCOL(X) <= 1){
    stop("X must have at least two columns.")
  }
  if(NROW(X) <= (NCOL(X) + intercept)){
    stop("X must have dimensions n > p (intercept=FALSE) or n > p+1 (intercept=TRUE)")
  }

  if(is.data.frame(X)){
    X.names <- names(X)
    X <- as.matrix(X, rownames.force = F)
  } else if(is.matrix(X)){
    X.names <- colnames(X)
  } else{
    stop("X must be a numeric matrix or data frame.")
  }

  n <- NROW(X)
  p <- NCOL(X)
  X.raw <- ifelse(n < 2*p+intercept, X, NA)

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
    if(n < 2*p+intercept){
      X <- rbind(X, matrix(0, 2*p+intercept - n, p))
    }
    kn_variables <- if("intercept" %in% names(formals(knockoffs))){
      knockoffs(X, intercept = intercept)
    } else{
      knockoffs(X)
    }
    X <- kn_variables$X
    X_kn <- kn_variables$Xk
  } else{
    stop("knockoff matrix/generating method of incorrect type.")
  }

  # record X and Xk before rotating
  X.org <- X
  X_kn.org <- X_kn

  # QR decomposition of [X, Xk]
  QR <- qr(cbind(X, X_kn))

  pivot_back <- sort.list(QR$pivot)

  Q_XXk <- qr.Q(QR, complete = F)
  R_XXk <- qr.R(QR, complete = F)[, pivot_back]

  # check regularity conditions
  tol <- max(abs(diag(R_XXk)[1:2*p])) * 1e-5
  if(max(abs(R_XXk[lower.tri(R_XXk, diag = F)])) > tol){
    stop("In the QR decompostion of the augmented maxtrix [X, Xk], ",
         "the resulted R matrix is not upper triangular. ",
         "Please contact the maintainer.")
  }
  if(min(abs(diag(R_XXk)[1:p])) < tol){
    stop(paste0("X", ifelse(intercept, "(along with intercept)", ""),
                " doesn't have full column rank."))
  }

  # rotate X and Xk by Q, then the new [X, Xk] are in the first 2p coordinates
  X <- R_XXk[, 1:p]
  X_kn <- R_XXk[, (p+1):(2*p)]

  # compute the a basic matrix needed by cknockoff
  # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
  vj_mat <- sapply(1:p, function(j){
    # Q_X[, j] = unit of X_j orthogonally projected onto X_{1:j-1}
    # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
    #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
    # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
    Q_X <- diag(2*p)[, 1:p]
    coord <- forwardsolve(t(R_XXk[j:p,j:p]), c(1, rep(0, p-j)))
    vj <- Q_X[, j:p] %*% matrix(coord, nrow = p-j+1)
    vj <- vj / sqrt(sum(vj^2))
  })

  X.pack.data <- list(X = X, X_kn = X_kn,
                      vj_mat = vj_mat,
                      X.raw = X.raw, X.org = X.org, X_kn.org = X_kn.org,
                      XXk.org.basis = Q_XXk)

  X.pack <- structure(c(X.pack.data, list(X.names = X.names, X.org.nrow = n)),
                      class = "cknockoff.X.pack")

  return(X.pack)

}



transform_y <- function(X.pack, y,
                        intercept,
                        randomize = F){
  n <- length(y)
  p <- NCOL(X.pack$X)

  # center y if necessary
  y.raw <- y
  y <- scale(y, center = intercept, scale = F)

  # compute RSS and degree of freedom in y ~ X
  RSS_X <- sum(y^2) - sum((matrix(y, nrow=1) %*% X.pack$XXk.org.basis[1:n, 1:p])^2)
  df_X <- n - p - intercept

  # augment y if needed
  if(n < 2*p+intercept){
    if(randomize){
      y.extra <- rnorm(2*p-n+intercept, sd = sqrt(RSS_X / df_X))
    } else{
      y.extra <- with_seed(0, rnorm(2*p-n+intercept, sd = sqrt(RSS_X / df_X)))
    }
    fitted_intercept <- if(intercept){
      unname(lm(y.raw ~ X.pack$X.raw + 1)$coefficients[1])
    } else{ 0 }

    y <- c(y.raw, y.extra + fitted_intercept)
    y <- scale(y, center = intercept, scale = F)
  }

  # record the original y before rotation
  y.org <- y

  # rotate y in the same way as we did for [X Xk]
  y.norm2 <- sum(y^2)
  y <- as.vector(matrix(y, nrow=1) %*% X.pack$XXk.org.basis)

  # the component of y not in the 2p-dim subspace (residue of y~[X Xk])
  # is recorded separately as RSS_XXk and df_XXk
  if(n < 2*p+intercept+1){
    RSS_XXk <- RSS_X
    df_XXk <- df_X
  } else{
    RSS_XXk <- y.norm2 - sum(y^2)
    df_XXk <- n - 2*p
  }

  y.data <- list(y = y, y.org = y.org,
                 RSS_X = RSS_X, RSS_XXk = RSS_XXk,
                 df_X = df_X, df_XXk = df_XXk)
  return(y.data)

}

# prepare the quantities needed by cknockoff based on y
process_y <- function(X.pack, y.data){
  p <- NCOL(X.pack$X)

  y_Pi_X <- c(y.data$y[1:p], rep(0, p))

  vjy_obs <- c(matrix(y.data$y, nrow=1) %*% X.pack$vj_mat)

  y_Pi_Xnoj <- sapply(1:p, function(j){
    y_Pi_X - X.pack$vj_mat[, j] * vjy_obs[j]
  })

  RSS_Xnoj <- y.data$RSS_X + vjy_obs^2
  # sigmahat.XXk_res <- sqrt(y.data$RSS_XXk / y.data$df_XXk)

  # placeholder, will be computed later in cknockoff()
  Xy_bias <- rep(NA, p)

  return(list(y.data = y.data, vjy_obs = vjy_obs,
              RSS_Xnoj = RSS_Xnoj,
              y_Pi_Xnoj = y_Pi_Xnoj,
              Xy_bias = Xy_bias))
}


# load packages and prepare snips for parallel computing
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

