## =====================================================================
##  The predictors of every data set (docs/validation_tests.md, "The data
##  sets"). Each data file draws them inside its own random stream
##  (R/shared/random_streams.R).
## =====================================================================


## n x k standard normal predictors with the same correlation between every
## pair, named x1 ... xk; centred if asked (the model conditions on them)
correlated_predictors <- function(n, k, correlation, centred = FALSE) {
  population <- matrix(correlation, k, k)
  diag(population) <- 1

  predictors <- matrix(rnorm(n * k), n) %*% chol(population)
  if (centred) predictors <- sweep(predictors, 2, colMeans(predictors))
  colnames(predictors) <- paste0("x", seq_len(k))
  predictors
}
