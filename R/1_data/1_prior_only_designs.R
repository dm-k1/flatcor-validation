## =====================================================================
##  The prior-only designs (study 1a): fixed predictors for drawing the
##  prior without data. Why these: docs/validation_tests.md, "The data sets".
##
##    one                1 predictor
##    three_correlated   3 predictors, correlation .4
##    eight              8 uncorrelated predictors
##  each n = 60, drawn once from its own stream.
## =====================================================================


data_prior_only_designs <- function(n = 60) {
  list(
    one = with_stream_of("prior-only-one", {
      correlated_predictors(n, 1, 0)
    }),
    three_correlated = with_stream_of("prior-only-three-correlated", {
      correlated_predictors(n, 3, 0.4)
    }),
    eight = with_stream_of("prior-only-eight", {
      correlated_predictors(n, 8, 0)
    })
  )
}

## The outcome the Stan model is given (it is ignored: prior_only = 1)
data_prior_only_outcome <- function(design_name, n = 60) {
  with_stream_of(paste0("prior-only-outcome-", design_name), rnorm(n))
}
