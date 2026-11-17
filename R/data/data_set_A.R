## =====================================================================
##  Data set A: fixed true effects (study 2b)
##
##  Predictors multivariate normal with unit variances and equal correlation
##  rho; noise SD 1, so slopes are in units of the noise SD. The slope of x1
##  is set by the design; of the other slopes, "sparse" has one at 0.3 and the
##  rest 0, "dense" has all at 0.1. Each dataset comes from its own random
##  stream, given by its ID (R/data/data_sets.R).
## =====================================================================

cell_id_A <- function(condition)
  with(condition, sprintf("A-k%d-n%d-rho%.1f-%s-slope%.1f", k, n, rho, others, slope_x1))

slopes_A <- function(k, others, slope_x1) {
  other_slopes <- if (others == "sparse") c(0.3, rep(0, k - 2)) else rep(0.1, k - 1)
  c(slope_x1, other_slopes)
}

## One dataset of set A: replication r of the cell in condition
generate_A <- function(condition, replication) {
  k <- condition$k; n <- condition$n
  covariance <- matrix(condition$rho, k, k); diag(covariance) <- 1
  with_stream_of(paste(cell_id_A(condition), replication, sep = "-r"), {
    predictors <- matrix(rnorm(n * k), n) %*% chol(covariance)
    colnames(predictors) <- paste0("x", seq_len(k))
    outcome <- drop(predictors %*% slopes_A(k, condition$others, condition$slope_x1)) + rnorm(n)
    list(predictors = predictors, outcome = outcome)
  })
}
