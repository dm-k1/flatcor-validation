## =====================================================================
##  Data for a smaller and a larger model (study 3b, items 3b.3-3b.6): the
##  rates of the Bayes factor as n grows. Why: docs/validation_tests.md,
##  "The data sets".
##
##  Standard normal predictors, noise SD 1. The larger model has the smaller
##  model's k_small predictors plus q more:
##    "x3 null next to two real"   slopes 0.3, 0.2 | 0         (the smaller model is true)
##    "x1 real"                    slope 0.5 | q x 0           (the smaller model is true)
##    "intercept only"             - | q x 0                   (the smaller model is true)
##    "alternative"                0.2, 0 | x1 = 0.15          (the larger model is true)
## =====================================================================


## One dataset: replication r of a cell (setting, q, n); returns the
## predictors of both models and the outcome
data_nested_models <- function(cell, replication) {
  k <- cell$k_small + cell$q
  id <- sprintf("nested-%s-q%d-n%d-r%d", cell$setting, cell$q, cell$n, replication)

  with_stream_of(id, {
    x <- correlated_predictors(cell$n, k, 0)
    slopes <- switch(cell$setting,
      `x3 null next to two real` = c(0.3, 0.2, 0),
      `intercept only` = rep(0, k),
      `x1 real` = c(0.5, rep(0, cell$q)),
      alternative = c(0.15, 0.2, 0)
    )
    smaller <- if (cell$setting == "alternative") {
      x[, 2:3]
    } else {
      x[, seq_len(cell$k_small), drop = FALSE]
    }
    list(smaller = smaller, larger = x, y = drop(x %*% slopes) + rnorm(cell$n))
  })
}
