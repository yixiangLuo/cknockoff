


y_to_t <- function(y, vj_mat, sigmahat){
  y <- matrix(y, nrow = 1)
  v_y <- y %*% vj_mat

  tvals <- c(v_y / sigmahat)

  return(tvals)
}

# compute the p-values for t-statistics
pvals_t <- function(tvals, df, side = "two"){
  if (side == "right"){
    pvals <- pt(tvals, df = df, lower.tail = FALSE)
  } else if (side == "left"){
    pvals <- pt(tvals, df = df, lower.tail = TRUE)
  } else if (side == "two"){
    pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
  }

  return(pvals)
}


# the complement of a union of disjoint intervals.
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
# univ_left, univ_right: the universal interval
interval_complement <- function(intervals, univ_left = 0, univ_right = 1){
  left <- c(univ_left, intervals$right)
  right <- c(intervals$left, univ_right)

  left_unique <- NULL
  right_unique <- NULL
  for(ind in 1:length(left)){
    if(left[ind] != right[ind]){
      left_unique <- c(left_unique, left[ind])
      right_unique <- c(right_unique, right[ind])
    }
  }

  complement <- list(left = left_unique, right = right_unique)

  return(complement)
}

# the intersection of two "union of disjoint intervals".
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_intersect <- function(intervals1, intervals2){
  left <- NULL
  right <- NULL
  intv1_num <- length(intervals1$left)
  intv2_num <- length(intervals2$left)
  if(min(intv1_num, intv2_num) > 0){
    for(ind1 in 1:intv1_num){
      for(ind2 in 1:intv2_num){
        low <- max(intervals1$left[ind1], intervals2$left[ind2])
        up <- min(intervals1$right[ind1], intervals2$right[ind2])
        if(low < up){
          left <- c(left, low)
          right <- c(right, up)
        }
      }
    }
  }
  intersect <- list(left = left, right = right)

  return(intersect)
}

# the union of two "union of disjoint intervals".
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_union <- function(intervals1, intervals2){
  complement1 <- interval_complement(intervals1, univ_left = -Inf, univ_right = Inf)
  complement2 <- interval_complement(intervals2, univ_left = -Inf, univ_right = Inf)

  complements_intersect <- interval_intersect(complement1, complement2)

  union <- interval_complement(complements_intersect, univ_left = -Inf, univ_right = Inf)

  return(union)
}

# intervals1 \setminus intervals2
# intervals$left: left bounds of the intervals in the increasing order;
# intervals$right: right bounds of the intervals in the increasing order
interval_minus <- function(intervals1, intervals2){
  if(is.null(intervals1$left)){
    return(list(left = NULL, right = NULL))
  } else if(is.null(intervals2$left)){
    return(list(left = intervals1$left, right = intervals1$right))
  } else{
    intv2_comp <- interval_complement(intervals2,
                                      univ_left = min(c(intervals1$left, intervals2$left)),
                                      univ_right = max(c(intervals1$right, intervals2$right)))
    minus <- interval_intersect(intervals1, intv2_comp)
    return(minus)
  }
}

# if x is in intervals
x_in_intervals <- function(x, intervals){
  if(is.null(intervals$left)){
    return(FALSE)
  }
  int_num <- length(intervals$left)
  for(i in 1:int_num){
    if(x >= intervals$left[i] & x <= intervals$right[i]){
      return(TRUE)
    }
  }
  return(FALSE)
}



# only take one iterator
iterate_seq <- function(...){
  arg <- list(...)[1]
  argname <- names(arg)
  arg <- arg[[1]]

  return(list(arg = arg, argname = argname))
}

`%do_seq%` <- function(obj, expr){
  result <- NULL
  if(length(obj$arg) > 0){
    for(item in obj$arg){
      envir <- parent.frame()
      envir[[obj$argname]] <- item
      result <- c(result, eval(substitute(expr), envir = envir))
    }
  }
  return(result)
}

# borrowed from knockoff
# Evaluate an expression with the given random seed, then restore the old seed
with_seed = function(seed, expr) {
  seed.old = if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed'))
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}

# borrowed from https://stackoverflow.com/a/43329945
match.call.defaults <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for(i in setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )


  match.call(sys.function(sys.parent()), call)
}
