# cknockoff

An R package for the cknockoff procedures


## Overview

The cKnockoff method is a multiple testing procedure applied to 
the linear model for variable selection with FDR control. It almost surely
dominates the fixed-X knockoff procedure.

This package implement the cKnockoff procedure efficiently.


## Installation

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("yixiangluo/cknockoff")
```

## Usage Examples

``` r
p <- 100; n <- 300; k <- 15
X <- matrix(rnorm(n*p), n)
nonzero <- sample(p, k)
beta <- 2.5 * (1:p %in% nonzero)
y <- X %*% beta + rnorm(n)
print(which(1:p %in% nonzero))

# Basic usage
result <- cknockoff(X, y, alpha = 0.05, n_cores = 1)
print(result$selected)

# knockoff rejection
library("knockoff")
kn.result <- knockoff.filter(X, y,
                             knockoffs = create.fixed,
                             statistic = stat.glmnet_coefdiff_lm,
                             fdr = 0.05)
print(kn.result$selected)

# improve knockoff result
result <- cknockoff(prelim_result = kn.result, n_cores = 2)
print(result$selected)

# improve previous cknockoff result
result <- cknockoff(prelim_result = result)
print(result$selected)

# improve previous cknockoff result with cknockoff*
result <- cknockoff(prelim_result = result, n_cores = 2, Rstar_refine = T)
print(result$selected)
```
