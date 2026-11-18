## =====================================================================
##  An ordinary regression dataset (study 1b), like the package's example. Why:
##  docs/validation_tests.md, "The data sets".
##
##    n = 60, three predictors with correlation .4, slopes 0.3, 0.2 and 0,
##    noise SD 1
## =====================================================================


data_ordinary_regression <- function(replication) {
  with_stream_of(sprintf("example-r%d", replication), {
    x <- correlated_predictors(60, 3, 0.4)
    y <- drop(x %*% c(0.3, 0.2, 0)) + rnorm(60)
    list(x = x, y = y)
  })
}
