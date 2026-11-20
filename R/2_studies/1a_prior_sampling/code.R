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
## percentile of every quantity has a Monte Carlo SE of at most Δ / 1.96,
## Δ = delta_share_1a x the quantity's SD (docs/validation_tests.md)
sample_setting_1a <- function(setting, predictors, model) {
  covariance <- covariance_of(predictors)
  stan_data <- list(
    N = nrow(predictors),
    K = ncol(predictors),
    X = predictors,
    y = data_prior_only_outcome(setting$design, nrow(predictors)),
    prior_only = 1L,
    proper_stand_ins = 1L,
    per_observation = 0L,   # no log_lik or y_rep: not needed here
    r_squared_shape = prior_shape(setting$prior, ncol(predictors))
  )
  name <- name_of("1a", setting$prior, setting$design)
  seed <- seed_of(name)
  started <- Sys.time()

  # the largest Monte Carlo SE of a percentile, as a share of its target
  # (the MCSE shrinks with the square root of the draws)
  share_of_target <- function(quantities) {
    max(vapply(quantities, function(x) {
      mcse <- vapply(
        percentiles_1a,
        function(p) posterior::mcse_quantile(x, probs = p),
        numeric(1)
      )
      max(mcse) / target_mcse(delta_share_1a * sd(as.vector(x)))
    }, numeric(1)))
  }

  # lengthen the chains until every percentile reaches its target
  iterations <- ceiling(start_draws_1a / 4)
  for (attempt in 1:6) {
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
    quantities <- quantities_1a(draws, covariance)
    smallest_ess <- min(vapply(quantities, posterior::ess_bulk, numeric(1)))
    share <- share_of_target(quantities)
    if (share <= 1) break
    iterations <- ceiling(iterations * 1.2 * share^2)
  }

  # the sampler per chain (CmdStan's diagnostics and times)
  diagnostics <- fit$diagnostic_summary(quiet = TRUE)
  times <- fit$time()$chains
  chains <- data.frame(
    chain = seq_along(diagnostics$num_divergent),
    divergent = diagnostics$num_divergent,
    max_treedepth_hits = diagnostics$num_max_treedepth,
    ebfmi = diagnostics$ebfmi,
    warmup_seconds = times$warmup,
    sampling_seconds = times$sampling
  )

  list(
    setting = setting,
    quantities = quantities_1a(draws, covariance),
    iterations = iterations,
    smallest_ess = smallest_ess,
    id = name,
    seed = seed,
    attempts = attempt,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    chains = chains
  )
}

## Layer 1 (R/shared/per_dataset.R): every prior draw of every quantity, one
## row per draw
draws_table_1a <- function(samples) {
  bind_rows(lapply(samples, function(sampled) {
    bind_rows(lapply(names(sampled$quantities), function(quantity) {
      draws <- sampled$quantities[[quantity]]
      tibble(
        prior = sampled$setting$prior,
        design = sampled$setting$design,
        quantity = quantity,
        chain = rep(seq_len(ncol(draws)), each = nrow(draws)),
        iteration = rep(seq_len(nrow(draws)), times = ncol(draws)),
        value = as.vector(draws)
      )
    }))
  }))
}

## Layer 1: each setting's run, one row per setting (its chains as a list
## column)
settings_table_1a <- function(samples) {
  bind_rows(lapply(samples, function(sampled) {
    tibble(
      prior = sampled$setting$prior,
      design = sampled$setting$design,
      id = sampled$id,
      seed = sampled$seed,
      iterations_per_chain = sampled$iterations,
      attempts = sampled$attempts,
      smallest_ess = sampled$smallest_ess,
      seconds = sampled$seconds,
      machine = machine_of(),
      chains = list(sampled$chains)
    )
  }))
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

## Layer 1: the recipe's draws for every prior (the three-predictor design),
## one row per draw
recipe_table_1a <- function(predictors) {
  covariance <- covariance_of(predictors)
  bind_rows(lapply(unique(settings_1a$prior), function(prior) {
    partials <- recipe_partials_1a(prior, covariance)
    bind_rows(lapply(seq_len(ncol(partials)), function(j) {
      tibble(
        prior = prior, design = "three_correlated",
        quantity = sprintf("partial correlation of x%d", j),
        draw = seq_len(nrow(partials)), value = partials[, j]
      )
    }))
  }))
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
  # Δ: 10% of the quantity's SD under the prior (the recipe's, or Beta(a, 1)'s)
  prior_sd <- if (startsWith(name, "partial")) {
    sd(reference)
  } else {
    sqrt(a / ((a + 1)^2 * (a + 2)))
  }

  quantiles <- expectation_rows(
    prior = setting$prior,
    design = setting$design,
    variable = name,
    percentile = 100 * percentiles_1a,
    quantity = paste("quantile of", name),
    expected = expected,
    estimate = estimate,
    mcse = mcse
  )
  quantiles$delta <- delta_share_1a * prior_sd

  list(
    quantiles = with_targets(quantiles),
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

## Every quantity of every setting, from the per-dataset files: the draws
## (`draws`, one row per draw), the recipe's draws (`recipe`) and the runs
## (`settings`); `predictors`: the designs' predictors
compare_1a <- function(draws, recipe, settings, predictors) {
  lapply(seq_len(nrow(settings)), function(i) {
    setting <- settings[i, c("prior", "design")]
    covariance <- covariance_of(predictors[[setting$design]])
    a <- prior_shape(setting$prior, ncol(covariance))
    of_setting <- draws[draws$prior == setting$prior & draws$design == setting$design, ]

    # the recipe's partial correlations as a draws x 3 matrix
    recipe_matrix <- if (setting$design == "three_correlated") {
      of_prior <- recipe[recipe$prior == setting$prior, ]
      sapply(sprintf("partial correlation of x%d", 1:3), function(q) {
        of_q <- of_prior[of_prior$quantity == q, ]
        of_q$value[order(of_q$draw)]
      })
    }

    rows <- lapply(unique(of_setting$quantity), function(name) {
      of_quantity <- of_setting[of_setting$quantity == name, ]
      x <- matrix(NA_real_, max(of_quantity$iteration), max(of_quantity$chain))
      x[cbind(of_quantity$iteration, of_quantity$chain)] <- of_quantity$value
      compare_quantity_1a(x, name, setting, a, recipe_matrix)
    })

    list(
      quantiles = do.call(rbind, lapply(rows, `[[`, "quantiles")),
      ks = do.call(rbind, lapply(rows, `[[`, "ks")),
      iterations = data.frame(
        prior = setting$prior,
        design = setting$design,
        iterations_per_chain = settings$iterations_per_chain[i],
        smallest_ess = round(settings$smallest_ess[i])
      )
    )
  })
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
