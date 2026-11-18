## =====================================================================
##  Data where the number of predictors grows with n (study 3b, items 3b.8
##  and 3b.9). Why: docs/validation_tests.md, "The data sets".
##
##  Two growth rates, k = n / 4 and k = sqrt(n), with k at most 90 (the range
##  in which bridge sampling has been shown to work). Standard normal
##  predictors, noise SD 1:
##    null true:        no predictor matters
##    small model true: 3 slopes of 0.3, the other k - 3 predictors null
## =====================================================================


## The cells: n for each growth rate
growing_k_cells <- function() {
  rbind(
    data.frame(n = c(50, 100, 200, 360), growth = "k = n / 4"),
    data.frame(n = c(50, 100, 200, 500, 1000, 2000), growth = "k = sqrt(n)")
  )
}

## k for a cell
k_of_growth <- function(n, growth) {
  if (growth == "k = n / 4") floor(n / 4) else floor(sqrt(n))
}

## One dataset: replication r of a cell, with truth "null" or "true"
data_growing_k <- function(cell, truth, replication) {
  k <- k_of_growth(cell$n, cell$growth)
  growth_name <- if (cell$growth == "k = n / 4") "quarter" else "root"
  id <- sprintf("growing-k-n%d-%s-%s-r%d", cell$n, growth_name, truth, replication)

  with_stream_of(id, {
    predictors <- correlated_predictors(cell$n, k, 0)
    slopes <- if (truth == "true") c(0.3, 0.3, 0.3, rep(0, k - 3)) else rep(0, k)
    list(predictors = predictors, outcome = drop(predictors %*% slopes) + rnorm(cell$n))
  })
}
