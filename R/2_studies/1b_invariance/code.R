## =====================================================================
##  Study 1b: the code
##  What is tested and why: docs/validation_tests.md, study 1b.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- The rescaled and recoded versions of a dataset (data: R/1_data/2_ordinary_regression.R) --

## Other units: x1 -> 100 x1 + 7, x2 -> -0.01 x2, x3 -> x3 - 50, y -> 1000 y - 3
rescaled_1b <- function(data) {
  list(
    x = cbind(100 * data$x[, 1] + 7, -0.01 * data$x[, 2], data$x[, 3] - 50),
    y = 1000 * data$y - 3
  )
}

## Other coding of the same model space: x1, x2 -> x1 + x2, x1 - x2
recoded_1b <- function(data) {
  list(
    x = cbind(data$x[, 1] + data$x[, 2], data$x[, 1] - data$x[, 2], data$x[, 3]),
    y = data$y
  )
}


## ---- The log Bayes factors of one dataset -----------------------------------------

## For one prior and one dataset: the five log Bayes factors, one row each
## with everything computed (R/shared/comparison_methods.R): the whole model
## against the intercept-only model on the original, rescaled and recoded
## data, and x1 (with it against without it) on the original and rescaled
## data. The changes are computed from these rows (changes_of_1b()).
invariance_1b <- function(prior, replication) {
  method <- names(contenders)[contenders == prior]
  original <- data_ordinary_regression(replication)
  versions <- list(
    original = original,
    rescaled = rescaled_1b(original),
    recoded = recoded_1b(original)
  )

  log_bf <- function(version, bayes_factor, without_columns, label) {
    data <- versions[[version]]
    started <- Sys.time()
    result <- log_bayes_factor_between(
      method,
      predictors_a = data$x,
      predictors_b = data$x[, without_columns, drop = FALSE],
      outcome = data$y,
      name = name_of("1b", prior, replication, label)
    )
    bind_cols(
      tibble(prior, replication, bayes_factor, version),
      select(result, -method),
      tibble(seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
             machine = machine_of())
    )
  }

  whole <- "whole model against intercept only"
  term <- "x1: with it against without it"
  bind_rows(
    # the whole model against the intercept-only model (no columns left)
    lapply(names(versions), function(v) log_bf(v, whole, 0, paste("whole", v))),
    # x1: the model with it against the model without it
    lapply(c("original", "rescaled"), function(v) log_bf(v, term, -1, paste("x1", v)))
  )
}

## One batch of datasets for one prior
invariance_batch_1b <- function(prior_row, replications) {
  bind_rows(lapply(replications, function(r) invariance_1b(prior_row$prior, r)))
}

## The change in log BF per dataset after rescaling and recoding, expected 0,
## from its rows: the changed log BF minus the original one. Its MCSE combines
## bridge sampling's reported errors of both (save_results_1b() corrects them
## by study 3a's part B).
changes_of_1b <- function(rows) {
  original <- rows[rows$version == "original", ]
  changed <- rows[rows$version != "original", ]
  paired <- merge(
    changed, original,
    by = c("prior", "replication", "bayes_factor"), suffixes = c("", "_original")
  )
  expectation_rows(
    prior = paired$prior,
    replication = paired$replication,
    bayes_factor = paired$bayes_factor,
    quantity = paste("change in log BF,", paired$version),
    expected = 0,
    estimate = paired$log_bf - paired$log_bf_original,
    mcse = sqrt(paired$error^2 + paired$error_original^2)
  ) |>
    arrange(prior, replication, bayes_factor, quantity)
}


## ---- Pooled over the datasets ------------------------------------------------------

## Per prior, Bayes factor and change: the mean change over the datasets,
## expected 0. Its MCSE is measured: the SD of the changes over the datasets
## / sqrt(n), which includes every source of noise, bridge sampling's among
## them (docs/validation_tests.md, "Monte Carlo errors are measured").
pooled_1b <- function(invariance) {
  invariance |>
    group_by(prior, bayes_factor, quantity) |>
    summarise(
      datasets = n(),
      expected = 0,
      mcse = sd(estimate) / sqrt(n()),
      estimate = mean(estimate),
      .groups = "drop"
    ) |>
    mutate(distance_in_mcse = (estimate - expected) / mcse) |>
    with_targets(deltas_1b)
}


## ---- The results file ---------------------------------------------------------

## `shortfall`: by how much bridge sampling's reported error falls short of
## the real spread (study 3a, part B). A change of one dataset is judged
## against the reported errors of its two Bayes factors times that factor;
## the reported errors are kept in column bridge_error.
save_results_1b <- function(invariance, shortfall) {
  pooled <- pooled_1b(invariance)
  invariance$bridge_error <- invariance$mcse
  invariance$mcse <- invariance$mcse * shortfall
  invariance$distance_in_mcse <- (invariance$estimate - invariance$expected) / invariance$mcse

  save_study_results_file("1b", list(
    pooled = pooled,
    invariance = invariance,
    beyond_reference = beyond_reference(invariance)
  ))
}
