## =====================================================================
##  Study 2b: the code
##  What is tested and why: docs/validation_tests.md, study 2b.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- One dataset -------------------------------------------------------------------

## Dataset `replication` of a cell (R/1_data/4_known_slopes.R), and for each method: x1's
## estimate and 95% interval, and the excess prediction error of the
## estimated slopes, (b - beta)' Sigma (b - beta)
estimates_2b <- function(cell, replication) {
  # the data, and the truth it came from
  data <- data_known_slopes(cell, replication)
  dataset_id <- paste(known_slopes_id(cell), replication, sep = "-r")
  true_slopes <- known_slopes(cell$k, cell$others, cell$slope_x1)
  covariance <- matrix(cell$rho, cell$k, cell$k)
  diag(covariance) <- 1

  columns <- colnames(data$predictors)
  frame <- data.frame(data$predictors, y = data$outcome)
  formula <- reformulate(columns, response = "y")

  # one row per method: x1's estimate and interval, and the prediction error
  summary_of <- function(method, slopes, lower, upper, divergent = NA) {
    error <- slopes - true_slopes
    data.frame(
      cell, replication, method,
      estimate = slopes[[1]],
      lower,
      upper,
      prediction_error = drop(t(error) %*% covariance %*% error),
      divergent
    )
  }
  from_draws <- function(method, draws, divergent = NA) {
    summary_of(
      method,
      slopes = colMeans(draws),
      lower = quantile(draws[, 1], 0.025, names = FALSE),
      upper = quantile(draws[, 1], 0.975, names = FALSE),
      divergent = divergent
    )
  }

  # the three priors: flatcor's Stan model
  priors <- lapply(names(contenders), function(method) {
    fit <- suppressWarnings(fit_regression(
      formula, frame,
      prior = contenders[[method]],
      chains = stan_draws_2b[["chains"]],
      iter_warmup = stan_draws_2b[["warmup"]],
      iter_sampling = stan_draws_2b[["sampling"]],
      seed = seed_of(dataset_id, method)
    ))
    draws <- as.matrix(as.data.frame(fit$draws)[, columns, drop = FALSE])
    from_draws(method, draws, fit$divergent)
  })

  # the references JZS and Zellner-Siow: BayesFactor's posterior sampler
  bayesfactor_draws <- function(scale, label) {
    with_stream_of(paste(dataset_id, label), {
      posterior <- BayesFactor::lmBF(
        formula, frame,
        rscaleCont = scale,
        posterior = TRUE,
        iterations = bayesfactor_draws_2b,
        progress = FALSE
      )
      as.matrix(as.data.frame(posterior)[, columns, drop = FALSE])
    })
  }
  references <- list(
    from_draws("JZS default", bayesfactor_draws(sqrt(2) / 4, "JZS")),
    from_draws("Zellner-Siow", bayesfactor_draws(1, "ZS"))
  )

  # ordinary least squares
  ols <- lm(formula, frame)
  ols_row <- summary_of(
    "OLS",
    slopes = coef(ols)[columns],
    lower = confint(ols)["x1", 1],
    upper = confint(ols)["x1", 2]
  )

  do.call(rbind, c(priors, references, list(ols_row)))
}

## One batch of datasets of one cell. `model`: flatcor's Stan model, compiled
## before, so that parallel workers find it instead of compiling it at once.
estimates_batch_2b <- function(cell, replications, model) {
  rows <- lapply(replications, function(r) estimates_2b(cell, r))
  do.call(rbind, rows)
}


## ---- Estimates against their references ----------------------------------------------

cell_label_2b <- function(rows) {
  sprintf("k = %d, n = %d, rho = %.1f, %s", rows$k, rows$n, rows$rho, rows$others)
}

## Bias of x1 (share of 0.2; reference 0) and coverage of its 95% interval
## (reference .95)
bias_and_coverage_2b <- function(estimates) {
  estimates |>
    mutate(cell = cell_label_2b(estimates)) |>
    group_by(k, n, rho, others, cell, method) |>
    summarise(
      datasets = n(),
      bias = (mean(estimate) - slope_x1[1]) / slope_x1[1],
      bias_mcse = sd(estimate) / sqrt(datasets) / slope_x1[1],
      coverage = mean(lower <= slope_x1[1] & slope_x1[1] <= upper),
      coverage_mcse = mcse_of_proportion(coverage, datasets),
      .groups = "drop"
    )
}

## Excess prediction error relative to OLS (reference 1), with a paired
## bootstrap MCSE (from its own fixed stream)
prediction_2b <- function(estimates) {
  ratios_in_cell <- function(rows, keys) {
    error_of <- function(name) {
      of_method <- rows[rows$method == name, ]
      of_method$prediction_error[order(of_method$replication)]
    }
    ols <- error_of("OLS")
    resamples <- replicate(
      if (smoke) 100 else 1000,
      sample(length(ols), replace = TRUE),
      simplify = FALSE
    )

    rows_per_method <- lapply(setdiff(methods_2b, "OLS"), function(name) {
      errors <- error_of(name)
      bootstrap <- vapply(resamples, function(r) mean(errors[r]) / mean(ols[r]), numeric(1))
      tibble(method = name, ratio = mean(errors) / mean(ols), ratio_mcse = sd(bootstrap))
    })
    bind_rows(rows_per_method)
  }

  with_stream_of("2b-bootstrap", {
    estimates |>
      mutate(cell = cell_label_2b(estimates)) |>
      group_by(k, n, rho, others, cell) |>
      group_modify(ratios_in_cell) |>
      ungroup()
  })
}

## One table: every estimate against its reference, with its Δ, target
## Monte Carlo SE and whether it is reached (R/shared/rounds.R)
results_table_2b <- function(bias_and_coverage, prediction) {
  bias <- transmute(
    bias_and_coverage, k, n, rho, others, cell, method,
    quantity = "bias of x1 (share of 0.2)",
    expected = 0, estimate = bias, mcse = bias_mcse
  )
  coverage <- transmute(
    bias_and_coverage, k, n, rho, others, cell, method,
    quantity = "coverage of x1's 95% interval",
    expected = 0.95, estimate = coverage, mcse = coverage_mcse
  )
  ratio <- transmute(
    prediction, k, n, rho, others, cell, method,
    quantity = "prediction error / OLS",
    expected = 1, estimate = ratio, mcse = ratio_mcse
  )

  bind_rows(bias, coverage, ratio) |>
    mutate(
      method = factor(method, levels = methods_2b),
      distance_in_mcse = (estimate - expected) / mcse
    ) |>
    with_targets(deltas_2b)
}

## Stan's divergent transitions, per prior
divergent_2b <- function(estimates) {
  estimates |>
    filter(method %in% names(contenders)) |>
    group_by(method) |>
    summarise(
      fits = n(),
      fits_with_divergences = sum(divergent > 0),
      divergent_transitions = sum(divergent)
    )
}


## ---- The results file ---------------------------------------------------------

save_results_2b <- function(estimates) {
  table <- results_table_2b(bias_and_coverage_2b(estimates), prediction_2b(estimates))

  save_study_results_file("2b", list(
    estimates = table,
    beyond_reference = beyond_reference(table, by = c("quantity", "method")),
    divergent = divergent_2b(estimates),
    per_dataset = estimates
  ))
}
