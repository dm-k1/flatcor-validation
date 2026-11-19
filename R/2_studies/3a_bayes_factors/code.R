## =====================================================================
##  Study 3a: the code
##  What is tested and why: docs/validation_tests.md, study 3a.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- One Bayes factor ---------------------------------------------------------

## Dataset `replication` of a cell, and the log Bayes factor for x1 under
## `prior`: the model with x1 against the model without it.
##   purpose  keeps the datasets of part A and part B apart
##   run      a repeat of the same dataset, with fresh posterior draws
bayes_factor_3a <- function(prior, k, n, replication, purpose, run = 1) {
  # the data: odd replications of part A come from the model without x1
  from_model_without_x1 <- purpose == "calibration" && replication %% 2 == 1
  data <- data_from_prior(
    prior, k, n,
    x1 = if (from_model_without_x1) "reduced" else "prior",
    purpose = purpose,
    replication = replication,
    location_scale = "fixed"
  )

  # the Bayes factor: each model fitted with flatcor's Stan model, the
  # marginal likelihoods by bridge sampling
  result <- log_bayes_factor_between(
    method = names(contenders)[contenders == prior],
    predictors_a = data$predictors,
    predictors_b = data$predictors[, -1, drop = FALSE],
    outcome = data$outcome,
    seed = seed_of("3a", purpose, prior, k, n, replication, run)
  )

  data.frame(
    prior, k, n, replication, run, from_model_without_x1,
    log_BF10 = result[["log_bf"]],
    error = result[["error"]],
    divergent = result[["divergent"]]
  )
}


## ---- Part A: calibration ------------------------------------------------------

## One batch of datasets of one cell
calibration_batch_3a <- function(cell, replications) {
  rows <- lapply(replications, function(r) {
    bayes_factor_3a(cell$prior, cell$k, cell$n, r, purpose = "calibration")
  })
  do.call(rbind, rows)
}

## Per bin of log BF10: the observed share of datasets without x1, against
## the mean of 1 / (1 + BF10). Its MCSE: the Monte Carlo SD of the share if
## the Bayes factors are calibrated.
calibration_bins_3a <- function(results) {
  bin_edges <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, Inf)

  # each dataset: its bin, and its predicted chance of having no x1
  binned <- results |>
    mutate(
      cell = sprintf("k = %d, n = %d", k, n),
      bin = cut(log_BF10, bin_edges),
      predicted = plogis(-log_BF10) # = 1 / (1 + BF10)
    )

  # per bin: observed share against the expected share
  binned |>
    group_by(prior, cell, bin) |>
    summarise(
      datasets = n(),
      observed = mean(from_model_without_x1),
      expected = mean(predicted),
      mcse = sqrt(sum(predicted * (1 - predicted))) / n(),
      .groups = "drop"
    ) |>
    filter(datasets >= bin_minimum_3a) |>
    transmute(
      prior, cell,
      bin = as.character(bin),
      datasets,
      quantity = "share without x1, per bin of log BF10",
      expected,
      estimate = observed,
      mcse,
      distance_in_mcse = (estimate - expected) / mcse
    )
}

## The data-averaged posterior, per prior and cell: the mean over the
## datasets of P(without x1 | data) = 1 / (1 + BF10), expected 1/2, the prior
## probability of the model without x1 (Schad et al., 2023). It decides part
## A: it uses every dataset's posterior probability in full. Its MCSE: the SD
## over the datasets / sqrt(datasets).
data_averaged_posterior_3a <- function(results) {
  results |>
    mutate(
      cell = sprintf("k = %d, n = %d", k, n),
      posterior_without_x1 = plogis(-log_BF10) # = 1 / (1 + BF10)
    ) |>
    group_by(prior, k, n, cell) |>
    summarise(
      datasets = n(),
      quantity = "data-averaged posterior",
      expected = 0.5,
      mcse = sd(posterior_without_x1) / sqrt(n()),
      estimate = mean(posterior_without_x1),
      .groups = "drop"
    ) |>
    mutate(distance_in_mcse = (estimate - expected) / mcse) |>
    with_targets(deltas_3a)
}

## The logistic calibration line of "without x1" on log BF01, per prior and
## cell: intercept expected 0, slope expected 1. It describes (no Δ): a line
## learns only from datasets whose Bayes factor is near 1, so with extreme
## Bayes factors it would need tens of thousands of datasets to reach one.
calibration_line_3a <- function(results) {
  fit_line <- function(rows, key) {
    # (huge Bayes factors give fitted probabilities of 0 or 1: harmless)
    fit <- suppressWarnings(glm(
      from_model_without_x1 ~ I(-log_BF10),
      family = binomial,
      data = rows
    ))
    coefficients <- summary(fit)$coefficients
    tibble(
      datasets = nrow(rows),
      quantity = c("calibration intercept", "calibration slope"),
      expected = c(0, 1),
      estimate = coefficients[, 1],
      mcse = coefficients[, 2]
    )
  }

  results |>
    mutate(cell = sprintf("k = %d, n = %d", k, n)) |>
    group_by(prior, k, n, cell) |>
    group_modify(fit_line) |>
    ungroup() |>
    mutate(distance_in_mcse = (estimate - expected) / mcse)
}


## ---- Part B: bridge sampling's error, from repeated runs ------------------------

## The same datasets, each computed `runs` times with fresh posterior draws
repeats_3a <- function(cell, datasets, runs) {
  rows <- lapply(datasets, function(r) {
    runs_of_r <- lapply(seq_len(runs), function(run) {
      bayes_factor_3a(cell$prior, cell$k, cell$n, r, purpose = "repeats", run = run)
    })
    do.call(rbind, runs_of_r)
  })
  do.call(rbind, rows)
}

## The spread over the runs (pooled SD) against the reported error (root mean
## square): expected ratio 1. Its MCSE: the relative SE of a pooled SD.
repeat_spread_3a <- function(repeats, runs) {
  # per dataset: the variance over its runs, and the mean squared reported error
  per_dataset <- repeats |>
    mutate(cell = sprintf("k = %d, n = %d", k, n)) |>
    group_by(prior, k, n, cell, replication) |>
    summarise(
      variance = var(log_BF10),
      mean_square_error = mean(error^2),
      .groups = "drop"
    )

  # per cell: pooled SD against the reported error
  per_dataset |>
    group_by(prior, k, n, cell) |>
    summarise(
      datasets = n(),
      pooled_sd = sqrt(mean(variance)),
      reported_error = sqrt(mean(mean_square_error)),
      .groups = "drop"
    ) |>
    mutate(
      quantity = "spread of repeated log BFs / reported error",
      expected = 1,
      estimate = pooled_sd / reported_error,
      mcse = estimate / sqrt(2 * datasets * (runs - 1)),
      distance_in_mcse = (estimate - expected) / mcse
    ) |>
    with_targets(deltas_3a)
}


## The factor by which bridge sampling's reported error falls short of the
## real spread, over every cell of part B: the pooled SD of the repeated runs
## against the root mean square reported error. Quantities about one given
## dataset (1b per dataset, 3b.1, 3b.2) multiply bridge sampling's error by
## it (docs/validation_tests.md, "Monte Carlo errors are measured").
bridge_shortfall_3a <- function(repeats) {
  per_dataset <- repeats |>
    group_by(prior, k, n, replication) |>
    summarise(
      variance = var(log_BF10),
      mean_square_error = mean(error^2),
      .groups = "drop"
    )
  sqrt(mean(per_dataset$variance)) / sqrt(mean(per_dataset$mean_square_error))
}


## ---- Also reported: bridge sampling's errors and Stan's divergent transitions ---

diagnostics_3a <- function(results) {
  results |>
    mutate(cell = sprintf("k = %d, n = %d", k, n)) |>
    group_by(prior, cell) |>
    summarise(
      datasets = n(),
      median_error = median(error),
      largest_error = max(error),
      fits_with_divergences = sum(divergent > 0),
      divergent_transitions = sum(divergent),
      .groups = "drop"
    )
}


## ---- The results file ---------------------------------------------------------

## Every table of the study in results/studies/study_3a.rds
## (R/3_figures/3a_bayes_factors.R makes the figures from it)
save_results_3a <- function(averaged, bins, line, spread, shortfall, diagnostics,
                            calibration, repeats) {
  distances <- rbind(
    averaged[, c("quantity", "distance_in_mcse")],
    bins[, c("quantity", "distance_in_mcse")],
    line[, c("quantity", "distance_in_mcse")],
    spread[, c("quantity", "distance_in_mcse")]
  )

  save_study_results_file("3a", list(
    data_averaged_posterior = averaged,
    calibration_bins = bins,
    calibration_line = line,
    repeated_runs = spread,
    bridge_shortfall = tibble(
      quantity = "real spread / reported error, over every cell of part B",
      factor = shortfall
    ),
    beyond_reference = beyond_reference(distances),
    diagnostics = diagnostics,
    bayes_factors = calibration,
    repeats = repeats
  ))
}
