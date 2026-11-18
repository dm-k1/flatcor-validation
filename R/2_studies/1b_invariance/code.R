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


## ---- The change in log Bayes factor -------------------------------------------------

## For one prior and one dataset: the change after rescaling and recoding,
## expected 0. Its MCSE combines bridge sampling's errors of both Bayes factors.
invariance_1b <- function(prior, replication) {
  method <- names(contenders)[contenders == prior]
  original <- data_ordinary_regression(replication)
  versions <- list(
    original = original,
    rescaled = rescaled_1b(original),
    recoded = recoded_1b(original)
  )

  log_bf <- function(version, without_columns, label) {
    data <- versions[[version]]
    log_bayes_factor_between(
      method,
      predictors_a = data$x,
      predictors_b = data$x[, without_columns, drop = FALSE],
      outcome = data$y,
      seed = seed_of("1b", prior, replication, label)
    )
  }

  # the whole model against the intercept-only model (no columns left)
  whole <- lapply(names(versions), function(v) log_bf(v, 0, paste("whole", v)))
  names(whole) <- names(versions)

  # x1: the model with it against the model without it
  term <- lapply(c("original", "rescaled"), function(v) log_bf(v, -1, paste("x1", v)))
  names(term) <- c("original", "rescaled")

  # the change, and its MCSE
  change <- function(pair, after) {
    c(
      estimate = pair[[after]][["log_bf"]] - pair$original[["log_bf"]],
      mcse = sqrt(pair[[after]][["error"]]^2 + pair$original[["error"]]^2)
    )
  }
  values <- rbind(
    change(whole, "rescaled"),
    change(whole, "recoded"),
    change(term, "rescaled")
  )

  expectation_rows(
    prior = prior,
    replication = replication,
    bayes_factor = c(
      "whole model against intercept only",
      "whole model against intercept only",
      "x1: with it against without it"
    ),
    quantity = paste("change in log BF,", c("rescaled", "recoded", "rescaled")),
    expected = 0,
    estimate = values[, "estimate"],
    mcse = values[, "mcse"]
  )
}

## One batch of datasets for one prior
invariance_batch_1b <- function(prior_row, replications) {
  rows <- lapply(replications, function(r) invariance_1b(prior_row$prior, r))
  do.call(rbind, rows)
}


## ---- Pooled over the datasets ------------------------------------------------------

## Per prior, Bayes factor and change: the mean change over the datasets,
## expected 0. Its MCSE combines the datasets' errors: sqrt(sum of MCSE^2) / n.
pooled_1b <- function(invariance) {
  invariance |>
    group_by(prior, bayes_factor, quantity) |>
    summarise(
      datasets = n(),
      expected = 0,
      estimate = mean(estimate),
      mcse = sqrt(sum(mcse^2)) / n(),
      .groups = "drop"
    ) |>
    mutate(distance_in_mcse = (estimate - expected) / mcse)
}


## ---- The results file ---------------------------------------------------------

save_results_1b <- function(invariance) {
  pooled <- pooled_1b(invariance)

  save_study_results_file("1b", list(
    pooled = pooled,
    invariance = invariance,
    beyond_reference = beyond_reference(invariance)
  ))
}
