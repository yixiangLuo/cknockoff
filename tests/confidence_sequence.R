
# test if we construct the confidence sequence correctly
# if so, the resulting "fails covering proportion" should be below alpha
test_confidence_sequence <- function(
  alpha, # significance level of the confidence sequence
  seq_length, # length of the data sequence
  n_test = 1000 # number of trials
  ){
  # number of trails that the confidence sequence doesn't cover the mean
  uncovered <- 0
  for(test_i in 1:n_test){
    # generate sequence of test data
    mean_vec <- c(-1, -1, 0)
    samples <- sample(mean_vec, size = seq_length, replace = T)
    samples <- samples + runif(seq_length, min = -1, max = 1)
    sample_bound <- list(lower = min(mean_vec)-1, upper = max(mean_vec)+1)
    true_mean <- mean(mean_vec)

    # linear transform the test sequence to [0, 1] using our function
    samples_unit <- sample_to_unit(samples, sample_bound)
    mean_unit <- sample_to_unit(true_mean, sample_bound)

    # call function m_in_HCS() to determine if the mean is covered
    mean_covered <- T
    for(step in 1:seq_length){
      if(!m_in_HCS(samples_unit[1:step], m = mean_unit, alpha)$in_HCS){
        mean_covered <- F
        break
      }
    }
    if(!mean_covered){
      uncovered <- uncovered + 1
    }
  }

  # print results
  print("-- test confidence sequence --")
  print(paste0("alpha: ", alpha,
               ", fails covering propotion: ", uncovered / n_test))
}


# test if our decision making function works properly
# if so
# 1. the "proportion of wrong decision when we have confidence" should be below
# "the effective significance level"
# 2. when "threshold" get farther from -2/3, the actual mean, the
# "proportion of decision with confidence" should be closer to 1
# 3. we expect to see a wrong decision w/o confidence rarely
test_decision <- function(
  rej_alpha, # significance level of rejecting the null
  accept_alpha, # significance level of accepting the null
  seq_length, # length of the data sequence, which has mean -2/3
  threshold, # we should reject the null if the mean -2/3 <= threshold
  n_test = 100 # number of trials
  ){
  # a matrix that record the rejecting/accepting decision
  # results[1,1]: number of accepting with confidence
  # results[1,2]: number of rejecting with confidence
  # results[2,1]: number of accepting without confidence
  # results[2,2]: number of rejecting without confidence
  results <- matrix(rep(0, 4), nrow = 2)
  for(test_i in 1:n_test){
    # generate sequence of test data
    mean_vec <- c(-1, -1, 0)
    samples <- sample(mean_vec, size = seq_length, replace = T)
    samples <- samples + runif(seq_length, min = -1, max = 1)
    sample_bound <- list(lower = min(mean_vec)-1, upper = max(mean_vec)+1)

    # initialize confidence to be FALSE
    confidence <- F
    for(step in 1:seq_length){
      # mimic the decision process in cknockoff,
      # as we see a sequence of f_j sample sequentially
      decision <- make_decision(samples[1:step], sample_bound, threshold, rej_alpha, accept_alpha)
      if(decision$confident){ # if we can make decision with confidence, do it
        results[1, as.integer(decision$reject)+1] <- 1 + results[1, as.integer(decision$reject)+1]
        confidence <- T
        break
      }
    }
    if(!confidence){ # if no confidence and we reach the end, make decision without confidence
      results[2, as.integer(decision$reject)+1] <- 1 + results[2, as.integer(decision$reject)+1]
    }
  }

  # compute the proportion
  results <- results / n_test
  if(mean(mean_vec) <= threshold){ # should we reject?
    alpha <- accept_alpha
    rej_type <- "accept"
    data <- results[, 1] / rowSums(results)
  } else{
    alpha <- rej_alpha
    rej_type <- "reject"
    data <- results[, 2] / rowSums(results)
  }

  print("-- test decision --")
  print(paste0("The wrong decision is: ", rej_type))
  print(paste0("Hence the effective significance level is: ", alpha))
  print(paste0("proportion of decision with confidence among all trials: ", sum(results[1, ])))
  print(paste0(" - proportion of wrong decision when we have confidence: ", data[1]))
  print(paste0("proportion of decision without confidence among all trials: ", sum(results[2, ])))
  print(paste0(" - proportion of wrong decision when we don't have confidence: ", data[2]))
}
