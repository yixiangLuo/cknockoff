# a technical utility, compute the denominator of lambdas in constructing the
# confidence sequence.
lambdas_denominator <- function(X_seq){
  t <- length(X_seq)
  if(t > 1){
    mu_hats <- (1/2 + cumsum(X_seq[1:(t-1)])) / (2:t)
    sigma2_hats <- (1/4 + cumsum((X_seq[1:(t-1)] - mu_hats)^2)) / (2:t)
    sigma2_hats <- c(1/4, sigma2_hats)
  } else{
    sigma2_hats <- 1/4
  }

  lambdas_deno <- sqrt((sigma2_hats * (1:t) * log((1:t) + 1)))

  return(lambdas_deno)
}

# check if m is in the Hedged Capital Confidence Sequence, see https://arxiv.org/abs/2010.09686
# sample mean is not guaranteed to be in the CI:
# set.seed(11); n=100000; x=runif(n); print(m_in_HCS(x, mean(x), alpha = 0.8))
# set.seed(11); n=100000; x=runif(n); print(m_in_HCS(x, 0.5, alpha = 0.8))
m_in_HCS <- function(X_seq, m, alpha = 0.05, c = 1/2, theta = 1/2, lambdas_deno = NULL){
  if(is.null(lambdas_deno)){
    lambdas_deno <- lambdas_denominator(X_seq)
  }

  lambdas <- sqrt(2*log(2/alpha)) / lambdas_deno

  lambdas_plus <- pmin(abs(lambdas), c/m)
  lambdas_minus <- pmin(abs(lambdas), c/(1-m))

  K_plus <- prod(1 + lambdas_plus * (X_seq - m))
  K_minus <- prod(1 - lambdas_minus * (X_seq - m))

  K_pm <- max(theta * K_plus, (1-theta) * K_minus)

  result <- list(in_HCS = (K_pm < 1 / alpha), K_val = K_pm)

  return(result)
}

# the lower and upper bound of the LHS of the inequality
get_Ej_bound <- function(alpha, p, weights, sample_coupling){
  if(!sample_coupling){
    upper <- max(1 * weights$rej_weight, 0 * weights$rest_weight)
    lower <- min((1/p - alpha) * weights$rej_weight,
                 (- alpha) * weights$rest_weight)
  } else{
    upper <- 1 * weights$rej_weight/2
    lower <- (1/p - alpha) * weights$rej_weight/2 -  alpha * weights$rest_weight/2
  }

  bounds <- list(upper = upper, lower = lower)

  return(bounds)
}

# linear transform the ineq = g_star - weight*alpha to [0,1]
sample_to_unit <- function(samples, sample_bound){
  width <- sample_bound$upper - sample_bound$lower
  samples <- (samples - sample_bound$lower) / width

  return(samples)
}

# decide whether to reject
decide_reject <- function(samples_unit, sample_mean, thr_unit, thr_Kval, alpha,
                          lambdas_deno, step_init = 1e-4, search_n = 10){
  # if the sample mean is in the HCCI, use it to make decision
  if(m_in_HCS(samples_unit, m = sample_mean, alpha = alpha , lambdas_deno = lambdas_deno)$in_HCS){
    return(sample_mean <= thr_unit)
  }
  # explore the local shape of the K function around thr_unit to decide whether the
  # confidence interval is on its left- or right-hand-side.
  # this is computationally more expensive and only used as a spare tire.
  else{
    step <- step_init
    for(i in 1:search_n){
      thr_left <- m_in_HCS(samples_unit, m = thr_unit-step, alpha = alpha , lambdas_deno = lambdas_deno)$K_val
      thr_right <- m_in_HCS(samples_unit, m = thr_unit+step, alpha = alpha , lambdas_deno = lambdas_deno)$K_val
      # confidence interval is on the LHS, reject
      if(thr_left < thr_Kval & thr_Kval < thr_right){
        return(T)
      }
      # confidence interval is on the RHS, don't reject
      if(thr_left > thr_Kval & thr_Kval > thr_right){
        return(F)
      }
      step <- step / 2
    }
    # don't reject (be conservative) if we cannot decide after running out of
    # computational budget
    return(F)
  }
}

# decide whether to reject and if we have confidence. Difference confidence level
# may be used for rejecting/not rejecting.
make_decision <- function(samples, sample_bound, threshold, rej_alpha, accept_alpha){
  # transform the samples and threshold to [0,1]
  samples_unit <- sample_to_unit(samples, sample_bound)
  thr_unit <- sample_to_unit(threshold, sample_bound)

  sample_mean <- mean(samples_unit)
  lambdas_deno <- lambdas_denominator(samples_unit)

  # use the sample mean if we can't decide with confidence
  confident <- F
  reject <- (sample_mean <= thr_unit)

  # make decision when we have confidence. We assume confidence level (1-rej_alpha)
  # for rejection is no smaller than the one (1-accept_alpha) for acception.

  # if the tighter confidence interval exclude the threshold
  thr_HCS_rej <- m_in_HCS(samples_unit, m = thr_unit, alpha = rej_alpha , lambdas_deno = lambdas_deno)
  if(!thr_HCS_rej$in_HCS){
    confident <- T
    reject <- decide_reject(samples_unit, sample_mean, thr_unit,
                            thr_Kval = thr_HCS_rej$K_val, alpha = rej_alpha, lambdas_deno = lambdas_deno)
  }
  else if(rej_alpha != accept_alpha){
    # if the looser confidence interval exclude the threshold (only decide with confidence
    # if it implies to accept)
    thr_HCS_acp <- m_in_HCS(samples_unit, m = thr_unit, alpha = accept_alpha , lambdas_deno = lambdas_deno)
    if(!thr_HCS_acp$in_HCS){
      if(!decide_reject(samples_unit, sample_mean, thr_unit,
                        thr_Kval = thr_HCS_acp$K_val, alpha = accept_alpha, lambdas_deno = lambdas_deno)){
        confident <- T
        reject <- F
      }
    }
  }

  decision <- list(reject = reject, confident = confident)

  return(decision)
}
