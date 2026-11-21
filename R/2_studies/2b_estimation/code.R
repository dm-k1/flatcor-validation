## =====================================================================
##  Study 2b: the code
##  What is tested and why: docs/validation_tests.md, study 2b.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- One dataset -------------------------------------------------------------------

## Dataset `replication` of a cell (R/1_data/4_known_slopes.R), and for each method:
## x1's estimate and 95% interval, the excess prediction error of the
## estimated slopes, (b - beta)' Sigma (b - beta), and everything else
## computed: per method the summary of every parameter (list column
## `variables`; for OLS its coefficient table), and for the three priors
## Stan's sampler per chain (`chains`), seeds and seconds
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
  sample_r_squared <- summary(lm(formula, frame))$r.squared

  # the summary of every parameter from posterior draws (the posterior package)
  summary_of_draws <- function(draws) {
    as.data.frame(posterior::summarise_draws(
      posterior::as_draws(draws),
      "mean", "median", "sd", "mad",
      ~ posterior::quantile2(.x, probs = c(0.025, 0.05, 0.95, 0.975)),
      "rhat", "ess_bulk", "ess_tail"
    ))
  }

  # one row per method: x1's estimate and interval, the prediction error, the rest
  row_of <- function(method, slopes, lower, upper, name, seed = NA, divergent = NA,
                     variables = NULL, chains = NULL, seconds = NA) {
    error <- slopes - true_slopes
    tibble(
      cell, replication, dataset = dataset_id, id = name, seed = seed, method,
      estimate = slopes[[1]],
      lower,
      upper,
      prediction_error = drop(t(error) %*% covariance %*% error),
      divergent,
      sample_r_squared = sample_r_squared,
      seconds = seconds,
      machine = machine_of(),
      variables = list(variables),
      chains = list(chains)
    )
  }
  from_draws <- function(method, draws, name, ...) {
    slopes <- as.matrix(as.data.frame(draws)[, columns, drop = FALSE])
    row_of(
      method,
      slopes = colMeans(slopes),
      lower = quantile(slopes[, 1], 0.025, names = FALSE),
      upper = quantile(slopes[, 1], 0.975, names = FALSE),
      name = name,
      variables = summary_of_draws(draws),
      ...
    )
  }

  # the three priors: flatcor's Stan model
  priors <- lapply(names(contenders), function(method) {
    name <- name_of(dataset_id, method)
    fit <- suppressWarnings(fit_regression(
      formula, frame,
      prior = contenders[[method]],
      chains = stan_draws_2b[["chains"]],
      iter_warmup = stan_draws_2b[["warmup"]],
      iter_sampling = stan_draws_2b[["sampling"]],
      seed = seed_of(name)
    ))
    from_draws(
      method, fit$draws, name,
      seed = fit$seed, divergent = fit$divergent, chains = fit$sampler,
      seconds = sum(unlist(fit$seconds))
    )
  })

  # the references JZS and Zellner-Siow: BayesFactor's posterior sampler
  bayesfactor_draws <- function(scale, label) {
    with_stream_of(paste(dataset_id, label), {
      started <- Sys.time()
      posterior <- BayesFactor::lmBF(
        formula, frame,
        rscaleCont = scale,
        posterior = TRUE,
        iterations = bayesfactor_draws_2b,
        progress = FALSE
      )
      draws <- as.matrix(posterior)
      attr(draws, "seconds") <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      draws
    })
  }
  references <- lapply(
    list(c("JZS default", "JZS", sqrt(2) / 4), c("Zellner-Siow", "ZS", 1)),
    function(reference) {
      draws <- bayesfactor_draws(as.numeric(reference[3]), reference[2])
      from_draws(
        reference[1], posterior::as_draws_matrix(draws),
        name = paste(dataset_id, reference[2]), seconds = attr(draws, "seconds")
      )
    }
  )

  # ordinary least squares, with its coefficient table
  ols <- lm(formula, frame)
  coefficients <- summary(ols)$coefficients
  ols_row <- row_of(
    "OLS",
    slopes = coef(ols)[columns],
    lower = confint(ols)["x1", 1],
    upper = confint(ols)["x1", 2],
    name = paste(dataset_id, "OLS"),
    variables = data.frame(
      variable = rownames(coefficients),
      mean = coefficients[, 1], sd = coefficients[, 2],
      t = coefficients[, 3], p = coefficients[, 4],
      q2.5 = confint(ols)[, 1], q97.5 = confint(ols)[, 2]
    )
  )

  bind_rows(c(priors, references, list(ols_row)))
}

## One batch of datasets of one cell. `model`: flatcor's Stan model, compiled
## before, so that parallel workers find it instead of compiling it at once.
estimates_batch_2b <- function(cell, replications, model) {
  compile_flatcor_models() # once per worker, one worker at a time (R/shared/pipeline.R)
  bind_rows(lapply(replications, function(r) estimates_2b(cell, r)))
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
## bootstrap MCSE. Each cell's bootstrap has its own random stream, named
## after the cell, so a cell's MCSE depends only on its own datasets, not on
## how many the other cells have (run_all.R judges each cell on its own, as
## soon as its batch is done)
prediction_2b <- function(estimates) {
  ratios_in_cell <- function(rows, keys) {
    error_of <- function(name) {
      of_method <- rows[rows$method == name, ]
      of_method$prediction_error[order(of_method$replication)]
    }
    ols <- error_of("OLS")
    resamples <- with_stream_of(paste("2b-bootstrap", keys$cell), {
      replicate(if (smoke) 100 else 1000, sample(length(ols), replace = TRUE), simplify = FALSE)
    })

    rows_per_method <- lapply(setdiff(methods_2b, "OLS"), function(name) {
      errors <- error_of(name)
      bootstrap <- vapply(resamples, function(r) mean(errors[r]) / mean(ols[r]), numeric(1))
      tibble(method = name, ratio = mean(errors) / mean(ols), ratio_mcse = sd(bootstrap))
    })
    bind_rows(rows_per_method)
  }

  estimates |>
    mutate(cell = cell_label_2b(estimates)) |>
    group_by(k, n, rho, others, cell) |>
    group_modify(ratios_in_cell) |>
    ungroup()
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

## `estimates`: the per-dataset rows, read from results/per_dataset/
## (study_2b_fits)
save_results_2b <- function(estimates) {
  table <- results_table_2b(bias_and_coverage_2b(estimates), prediction_2b(estimates))

  save_study_results_file("2b", list(
    estimates = table,
    beyond_reference = beyond_reference(table, by = c("quantity", "method")),
    divergent = divergent_2b(estimates)
  ))
}
