## =====================================================================
##  Data with known slopes (study 2b): fixed true effects. Why, and every
##  choice: docs/validation_tests.md, "The data sets".
##
##  Predictors with unit variances and the same correlation rho; noise SD 1,
##  so slopes are in units of the noise SD. The slope of x1 is set by the
##  cell; the other slopes are "sparse" (one at 0.3, the rest 0) or "dense"
##  (all 0.1).
## =====================================================================


## The true slopes of a cell
known_slopes <- function(k, others, slope_x1) {
  other_slopes <- if (others == "sparse") c(0.3, rep(0, k - 2)) else rep(0.1, k - 1)
  c(slope_x1, other_slopes)
}

## The ID of a cell (part of every dataset's ID)
known_slopes_id <- function(cell) {
  with(cell, sprintf("known-slopes-k%d-n%d-rho%.1f-%s-slope%.1f", k, n, rho, others, slope_x1))
}

## One dataset: replication r of a cell
data_known_slopes <- function(cell, replication) {
  with_stream_of(paste(known_slopes_id(cell), replication, sep = "-r"), {
    predictors <- correlated_predictors(cell$n, cell$k, cell$rho)
    slopes <- known_slopes(cell$k, cell$others, cell$slope_x1)
    outcome <- drop(predictors %*% slopes) + rnorm(cell$n)
    list(predictors = predictors, outcome = outcome)
  })
}
