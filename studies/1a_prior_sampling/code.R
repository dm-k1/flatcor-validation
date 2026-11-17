## =====================================================================
##  Study 1a: the code
##  What is tested and why: docs/validation_tests.md, study 1a.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- Quantities from the draws -------------------------------------------------

## The predictors' covariance (as the Stan model computes it)
covariance_of <- function(predictors) {
  crossprod(scale(predictors, scale = FALSE)) / nrow(predictors)
}

## The draws of one coefficient, or of sigma, as an iterations x chains matrix
coefficient_draws <- function(draws, j) {
  as.matrix(posterior::extract_variable_array(draws, sprintf("b[%d]", j))[, , 1])
}
sigma_draws <- function(draws) {
  as.matrix(posterior::extract_variable_array(draws, "sigma")[, , 1])
}

## The quantities whose distribution is checked, each an iterations x chains
## matrix: r^2 (one predictor) or R-squared, and with three predictors each
## coefficient's partial correlation
quantities_1a <- function(draws, covariance) {
  k <- ncol(covariance)
  chains <- ncol(sigma_draws(draws))

  # R-squared = f^2 / (1 + f^2), with f^2 = b' Sigma b / sigma^2
  coefficients <- sapply(seq_len(k), function(j) as.vector(coefficient_draws(draws, j)))
  coefficients <- matrix(coefficients, ncol = k)
  f_squared <- rowSums((coefficients %*% covariance) * coefficients) /
    as.vector(sigma_draws(draws))^2
  r_squared <- matrix(f_squared / (1 + f_squared), ncol = chains)

  if (k == 1) {
    return(list(`r^2` = r_squared))
  }

  # partial correlation: the coefficient in noise units times the predictor's
  # spread left after the others, rho = f_j / sqrt(1 + f_j^2)
  partials <- NULL
  if (k == 3) {
    spread_left <- sqrt(1 / diag(solve(covariance)))
    partials <- lapply(1:3, function(j) {
      f_j <- coefficient_draws(draws, j) / sigma_draws(draws) * spread_left[j]
      f_j / sqrt(1 + f_j^2)
    })
    names(partials) <- sprintf("partial correlation of x%d", 1:3)
  }

  c(list(`R-squared` = r_squared), partials)
}


## ---- The prior draws of one setting ------------------------------------------------

## flatcor's Stan model with the data switched off (prior_only = 1, with
## proper stand-ins for the intercept and the noise), lengthened until every
## quantity reaches the target effective sample size
sample_setting_1a <- function(setting, predictors, model) {
  covariance <- covariance_of(predictors)
  stan_data <- list(
    N = nrow(predictors),
    K = ncol(predictors),
    X = predictors,
    y = data_prior_only_outcome(setting$design, nrow(predictors)),
    prior_only = 1L,
    proper_stand_ins = 1L,
    r_squared_shape = prior_shape(setting$prior, ncol(predictors))
  )
  seed <- seed_of("1a", setting$prior, setting$design)

  # lengthen the chains until the smallest ESS reaches the target
  iterations <- ceiling(1.2 * target_ess_1a / 4)
  for (attempt in 1:4) {
    fit <- model$sample(
      data = stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = iterations,
      seed = seed,
      refresh = 0,
      show_messages = FALSE
    )
    draws <- fit$draws(c("b", "sigma"))
    ess <- vapply(quantities_1a(draws, covariance), posterior::ess_bulk, numeric(1))
    smallest_ess <- min(ess)
    if (smallest_ess >= target_ess_1a) break
    iterations <- ceiling(iterations * 1.2 * target_ess_1a / smallest_ess)
  }

  list(
    setting = setting,
    quantities = quantities_1a(draws, covariance),
    iterations = iterations,
    smallest_ess = smallest_ess
  )
}


## ---- Quantiles against their expected values ---------------------------------------

## The prior as a recipe (draw_coefficients_from_prior()), as partial
## correlations: the second route for the three-predictor design
recipe_partials_1a <- function(prior, covariance) {
  draws <- with_stream_of(paste0("prior-recipe-", prior, "-three"), {
    draw_coefficients_from_prior(reference_size_1a, covariance, prior = prior)
  })
  f <- sweep(draws, 2, sqrt(1 / diag(solve(covariance))), "*")
  f / sqrt(1 + f^2)
}

## One quantity: its quantiles against the expected ones, and its KS distance
## to the expected distribution next to the .01 critical value for its ESS
compare_quantity_1a <- function(x, name, setting, a, recipe) {
  estimate <- quantile(as.vector(x), percentiles_1a, names = FALSE)
  mcse_of <- function(values) {
    vapply(percentiles_1a, function(p) posterior::mcse_quantile(values, probs = p), numeric(1))
  }
  mcse <- mcse_of(x)
  ess <- posterior::ess_bulk(x)

  if (startsWith(name, "partial")) {
    # against the recipe, whose own MCSE is added
    reference <- recipe[, as.integer(sub(".*x", "", name))]
    expected <- quantile(reference, percentiles_1a, names = FALSE)
    mcse <- sqrt(mcse^2 + mcse_of(reference)^2)
    ks_distance <- suppressWarnings(ks.test(as.vector(x), reference)$statistic)
    critical <- 1.63 * sqrt((ess + reference_size_1a) / (ess * reference_size_1a))
    against <- "the recipe"
  } else {
    # against Beta(a, 1)
    expected <- qbeta(percentiles_1a, a, 1)
    ks_distance <- suppressWarnings(ks.test(as.vector(x), function(v) pbeta(v, a, 1))$statistic)
    critical <- 1.63 / sqrt(ess)
    against <- sprintf("Beta(%g, 1)", a)
  }

  list(
    quantiles = expectation_rows(
      prior = setting$prior,
      design = setting$design,
      variable = name,
      percentile = 100 * percentiles_1a,
      quantity = paste("quantile of", name),
      expected = expected,
      estimate = estimate,
      mcse = mcse
    ),
    ks = data.frame(
      prior = setting$prior,
      design = setting$design,
      variable = name,
      against = against,
      ks_distance = unname(ks_distance),
      critical_at_01 = critical,
      ess = round(ess)
    )
  )
}

## Every quantity of one setting
compare_setting_1a <- function(sampled, predictors) {
  setting <- sampled$setting
  covariance <- covariance_of(predictors)
  a <- prior_shape(setting$prior, ncol(covariance))
  recipe <- if (setting$design == "three_correlated") {
    recipe_partials_1a(setting$prior, covariance)
  }

  rows <- lapply(names(sampled$quantities), function(name) {
    compare_quantity_1a(sampled$quantities[[name]], name, setting, a, recipe)
  })

  list(
    quantiles = do.call(rbind, lapply(rows, `[[`, "quantiles")),
    ks = do.call(rbind, lapply(rows, `[[`, "ks")),
    iterations = data.frame(
      prior = setting$prior,
      design = setting$design,
      iterations_per_chain = sampled$iterations,
      smallest_ess = round(sampled$smallest_ess)
    )
  )
}


## ---- The results file ---------------------------------------------------------

save_results_1a <- function(comparisons) {
  quantiles <- do.call(rbind, lapply(comparisons, `[[`, "quantiles"))

  save_study_results_file("1a", list(
    quantiles = quantiles,
    beyond_reference = beyond_reference(quantiles),
    ks = do.call(rbind, lapply(comparisons, `[[`, "ks")),
    iterations = do.call(rbind, lapply(comparisons, `[[`, "iterations"))
  ))
}
