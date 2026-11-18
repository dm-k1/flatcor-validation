## =====================================================================
##  Data from the prior (studies 2a and 3a): the effect drawn from the prior
##  being tested. Why, and every choice: docs/validation_tests.md, "The data
##  sets".
##
##  One dataset = predictors x coefficients (from the prior) + noise
##
##    shared by the three priors (common random numbers): the predictors
##      (correlation .3, centred), the noise draws, and the intercept and
##      noise SD (drawn from the Stan model's proper stand-ins, N(0, 10) and
##      lognormal(0, 1), or fixed at 0 and 1)
##    from each prior: the coefficients, R-squared ~ Beta(a, 1) and a uniform
##      direction; with x1 = "reduced", x1 = 0 and the others from the prior
##      of the model without x1
## =====================================================================


## The coefficients in noise units: all from the prior, or x1 = 0 and the
## others from the prior of the model without x1
coefficients_from_prior <- function(prior, covariance, x1) {
  if (x1 == "prior") {
    return(drop(draw_coefficients_from_prior(1, covariance, prior = prior)))
  }
  k <- ncol(covariance)
  if (k == 1) {
    return(0)
  }
  others <- covariance[-1, -1, drop = FALSE]
  c(0, drop(draw_coefficients_from_prior(1, others, prior = prior)))
}

## One dataset. `purpose` (the study's use) is part of the ID, so the uses do
## not share datasets by accident. `location_scale`: "drawn" (2a, where sigma
## is a parameter checked) or "fixed" (3a, where the Bayes factor does not
## depend on it).
data_from_prior <- function(prior, k, n, x1 = "prior", purpose, replication,
                            location_scale = "drawn", predictor_correlation = 0.3) {
  # shared by the three priors
  shared_id <- sprintf("from-prior-%s-shared-k%d-n%d-r%d", purpose, k, n, replication)
  shared <- with_stream_of(shared_id, {
    list(
      predictors = correlated_predictors(n, k, predictor_correlation, centred = TRUE),
      noise = rnorm(n),
      intercept = rnorm(1, 0, 10),
      noise_sd = rlnorm(1, 0, 1)
    )
  })
  if (location_scale == "fixed") {
    shared$intercept <- 0
    shared$noise_sd <- 1
  }
  covariance <- crossprod(shared$predictors) / n

  # from each prior: the coefficients (noise units, then data units)
  prior_id <- sprintf("from-prior-%s-%s-k%d-n%d-x1_%s-r%d", purpose, prior, k, n, x1, replication)
  in_noise_units <- with_stream_of(prior_id, coefficients_from_prior(prior, covariance, x1))
  coefficients <- shared$noise_sd * in_noise_units

  outcome <- shared$intercept +
    drop(shared$predictors %*% coefficients) +
    shared$noise_sd * shared$noise
  f_squared <- sum((in_noise_units %*% covariance) * in_noise_units)

  list(
    predictors = shared$predictors,
    outcome = outcome,
    intercept = shared$intercept,
    noise_sd = shared$noise_sd,
    coefficients = coefficients,
    r_squared = f_squared / (1 + f_squared)
  )
}
