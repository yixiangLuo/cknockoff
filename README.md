# cknockoff

An R package for the cknockoff procedures


## Overview



## Installation

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("yixiangluo/cknockoff")
```

## Usage Examples

For details type `?cknockoff`.

``` r
# generate example data
p <- 100; n <- 300; k <- 15
X <- matrix(rnorm(n*p), n)
nonzero <- sample(p, k)
beta <- 3.5 * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)


# Basic usage
library("cknockoff")
result <- cknockoff(X, y, alpha = 0.05, n_cores = 1)
print(result$selected)


# improve knockoff and previous cknockoff result
library("knockoff")
kn.result <- knockoff.filter(X, y,
                             knockoffs = create.fixed,
                             statistic = stat.glmnet_coefdiff_lm,
                             fdr = 0.05)
print(kn.result$selected)

result <- cknockoff(prelim_result = kn.result, n_cores = 2)
print(result$selected)

result <- cknockoff(prelim_result = result)
print(result$selected)
```
