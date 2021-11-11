# weighted BH method
BH_weighted <- function(pvals, alpha,
                        weights = rep(1, length(pvals)) / length(pvals)){
  n <- length(pvals)

  adjust_pvals <- sort(pvals / weights) / (1:n)
  nrejs <- max(0, which(adjust_pvals <= alpha))

  rejs <- which(pvals <= nrejs * alpha * weights)

  return(list(nrejs = nrejs, rejs = rejs))
}

# q-values of BH
qvals_BH <- function(pvals){
  p <- length(pvals)
  fac <- 1:p
  ord <- order(pvals)

  adjust_pvals <- cummin(rev(pvals[ord] * p / fac))
  qvals <- rep(NA, p)
  qvals[ord] <- rev(adjust_pvals)

  return(qvals)
}

# pick the promising features from the unchecked ones
cKn_candidates <- function(kn.stats, pvals, alpha, record, selected){
  # when cknockoff is run the first time
  if(is.null(record) || record$iteration == 0) {
    # pick those selected by BH and having small p-values
    candidates <- intersect(which(pvals <= alpha/2), BH_weighted(pvals, 4*alpha)$rejs)
    # pick those having large knockoff W-statistics (in absolute value)
    kn_cand_num <- max(1, length(candidates), sum(abs(kn.stats) >=  kn.select(kn.stats, 1.5*alpha)$W_k_hat))
    candidates <- union(candidates, order(abs(kn.stats), decreasing = T)[1:kn_cand_num])
    candidates <- setdiff(candidates, selected)
  }
  # when we are recursively improving previous cknockoff result
  else {
    # rank the features based on similar creterion as above
    rank_p <- rank(pmax(pvals * 2, qvals_BH(pvals) / 4))
    rank_kn <- rank(-abs(kn.stats))
    ranks <- rank(rank_p + rank_kn, ties.method = "first")
    candidates <- 1:length(pvals)
    candidates[ranks] <- candidates

    # weed out those already checked
    candidates <- setdiff(candidates, record$checked_so_far)
    candidates <- setdiff(candidates, selected)
    # pick the top ones from the ranked list
    next_check <- min(length(candidates), record$next_check_num)
    if(next_check > 0){
      candidates <- candidates[1:next_check]
    }
  }

  # candidates <- sapply(candidates, function(candidate){
  #   prefer <- kn_prefer(candidate, kn.stats, alpha)
  #   if(prefer) return(candidate)
  #   else return(NA)
  # })
  # candidates <- candidates[!is.na(candidates)]

  return(candidates)
}
