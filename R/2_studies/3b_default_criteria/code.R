## =====================================================================
##  Study 3b: the code
##  What is tested and why: docs/validation_tests.md, study 3b.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- Small helpers -----------------------------------------------------------------

## The method's name for a prior ("flat_r" -> "flat R")
method_of <- function(prior) names(contenders)[contenders == prior]

## One Bayes factor as the row of its dataset: the design columns `design`
## (a one-row data frame), then everything computed (R/shared/comparison_methods.R),
## then the seconds it took and the machine
bayes_factor_row <- function(design, method, predictors_a, predictors_b, outcome, name) {
  started <- Sys.time()
  result <- log_bayes_factor_between(method, predictors_a, predictors_b, outcome, name)
  bind_cols(
    as_tibble(design),
    select(result, -any_of(names(design))),
    tibble(seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
           machine = machine_of())
  )
}

## One batch of datasets of one cell; `one` computes one dataset
batch_3b <- function(cell, replications, one) {
  bind_rows(lapply(replications, function(r) one(cell, r)))
}

## Slope of the median log BF against log n, with a bootstrap MCSE from
## resampling the datasets within each n
slope_of_medians <- function(rows) {
  by_n <- split(rows$log_bf, rows$n)
  sizes <- as.numeric(names(by_n))
  slope <- function(medians) coef(lm(medians ~ log(sizes)))[[2]]

  bootstrap <- replicate(bootstrap_replicates_3b, {
    resampled_medians <- vapply(by_n, function(v) median(sample(v, replace = TRUE)), numeric(1))
    slope(resampled_medians)
  })

  tibble(estimate = slope(vapply(by_n, median, numeric(1))), mcse = sd(bootstrap))
}


## ---- 3b.1, 3b.2: datasets with a set sample R-squared (R/1_data/5_set_r_squared.R) ---

## The log BF of the model against the intercept-only model, for one row of
## the design
formula_check_3b <- function(check) {
  data <- data_set_r_squared(check$n, check$k, check$r_squared)
  bayes_factor_row(
    check, method_of(check$prior),
    predictors_a = data$x,
    predictors_b = data$x[, 0, drop = FALSE],
    outcome = data$y,
    name = name_of("3b-formula", check$prior, check$k, check$r_squared)
  )
}


## ---- 3b.3-3b.6: the rates ----------------------------------------------------------

## Dataset `replication` of a rate cell (R/1_data/6_nested_models.R), and the
## log BF of the smaller model against the larger one (log_bf_smaller) and
## for the true model (log_bf): the truth is the smaller model, except under
## the alternative
rate_bayes_factor_3b <- function(cell, replication) {
  data <- data_nested_models(cell, replication)
  row <- bayes_factor_row(
    cbind(cell, replication = replication), method_of(cell$prior),
    predictors_a = data$smaller,
    predictors_b = data$larger,
    outcome = data$y,
    name = name_of("3b-rate", cell$setting, cell$q, cell$n, replication, cell$prior)
  )
  row$truth <- if (cell$setting == "alternative") "larger model" else "smaller model"
  row$log_bf_smaller <- row$log_bf
  if (cell$setting == "alternative") row$log_bf <- -row$log_bf
  row
}

## 3b.3-3b.5: the slopes against their asymptotic values, q/2 (or a against
## the intercept-only model)
rate_slopes_3b <- function(rates) {
  with_stream_of("3b-bootstrap-rates", {
    rates |>
      filter(setting != "alternative") |>
      group_by(item, setting, q, prior) |>
      group_modify(function(rows, keys) slope_of_medians(rows)) |>
      ungroup() |>
      mutate(
        quantity = "slope of median log BF against log n",
        expected = ifelse(setting == "intercept only", mapply(prior_shape, prior, q), q / 2),
        distance_in_mcse = (estimate - expected) / mcse,
        delta = relative_delta_3b * expected
      ) |>
      with_targets()
  })
}

## 3b.6: the slope of the median log BF10 against log n for a real effect
## (expected positive, a direction: its Δ by the relative rule, direction_delta())
alternative_slopes_3b <- function(rates) {
  alternative <- filter(rates, setting == "alternative")
  # none in the smoke test, whose corner leaves the alternative out (it
  # starts at n = 200): an empty table, as dplyr would otherwise fit a line
  # through no points
  if (nrow(alternative) == 0) {
    return(tibble(
      item = character(0), setting = character(0), q = numeric(0), prior = character(0),
      estimate = numeric(0), mcse = numeric(0), quantity = character(0),
      expected = character(0), delta = numeric(0), target_mcse = numeric(0),
      reached = logical(0)
    ))
  }
  with_stream_of("3b-bootstrap-alternative-slope", {
    alternative |>
      group_by(item, setting, q, prior) |>
      group_modify(function(rows, keys) slope_of_medians(rows)) |>
      ungroup() |>
      mutate(
        quantity = "slope of median log BF10 against log n",
        expected = "positive",
        delta = direction_delta(estimate, relative_delta_3b, direction_delta_3b)
      ) |>
      with_targets()
  })
}

## 3b.6: the median log BF10 per n (expected to increase with n)
alternative_3b <- function(rates) {
  bootstrap_median <- function(values) {
    sd(replicate(bootstrap_replicates_3b, median(sample(values, replace = TRUE))))
  }

  with_stream_of("3b-bootstrap-alternative", {
    rates |>
      filter(setting == "alternative") |>
      group_by(prior, n) |>
      summarise(
        median_log_bf = median(log_bf),
        mcse = bootstrap_median(log_bf),
        .groups = "drop"
      ) |>
      mutate(expected = "increases with n")
  })
}


## ---- 3b.7: intrinsic consistency -----------------------------------------------------

## The draws of x1's coefficient from the prior with the population
## covariance (three predictors, correlation .5; the limit) and with a sample
## covariance at each n: one row per draw (layer 1)
intrinsic_draws_of_3b <- function() {
  population <- matrix(0.5, 3, 3)
  diag(population) <- 1

  per_prior <- lapply(priors_3b, function(p) {
    limit <- with_stream_of(paste0("I5-", p, "-limit"), {
      draw_coefficients_from_prior(intrinsic_draws_3b, population, prior = p)[, 1]
    })
    per_n <- lapply(intrinsic_n_3b, function(n) {
      with_stream_of(sprintf("I5-%s-n%d", p, n), {
        x <- matrix(rnorm(n * 3), n) %*% chol(population)
        sample_covariance <- crossprod(scale(x, scale = FALSE)) / n
        draws <- draw_coefficients_from_prior(intrinsic_draws_3b, sample_covariance, prior = p)[, 1]
        tibble(prior = p, n = n, covariance = "sample", draw = seq_along(draws), value = draws)
      })
    })
    bind_rows(
      tibble(prior = p, n = NA_real_, covariance = "population", draw = seq_along(limit),
             value = limit),
      per_n
    )
  })
  bind_rows(per_prior)
}

## The prior with the sample covariance against that with the population
## covariance, per prior and n, by a two-sample KS distance (from the draws'
## file)
intrinsic_3b <- function(draws) {
  draws |>
    filter(covariance == "sample") |>
    group_by(prior, n) |>
    group_modify(function(rows, keys) {
      limit <- draws$value[draws$covariance == "population" & draws$prior == keys$prior]
      tibble(
        ks_distance = suppressWarnings(unname(ks.test(rows$value, limit)$statistic)),
        critical_at_01 = 1.63 * sqrt(2 / length(limit))
      )
    }) |>
    ungroup() |>
    mutate(expected = "falls below the critical value as n grows")
}


## ---- 3b.8-3b.9: k growing with n -----------------------------------------------------

## For the truth (R/1_data/7_growing_k.R): the null against the k-predictor
## model, or the true 3-predictor model against the full one
growing_bayes_factor_3b <- function(cell, replication) {
  data <- data_growing_k(cell, cell$truth, replication)
  smaller <- if (cell$truth == "null") {
    data$predictors[, 0, drop = FALSE]
  } else {
    data$predictors[, 1:3]
  }

  bayes_factor_row(
    cbind(cell, replication = replication), cell$method,
    predictors_a = smaller,
    predictors_b = data$predictors,
    outcome = data$outcome,
    name = name_of("3b-growing", cell$growth, cell$n, cell$truth, replication, cell$method)
  )
}

## The slopes, expected positive (a direction: its Δ by the relative rule,
## direction_delta())
growing_slopes_3b <- function(growing) {
  with_stream_of("3b-bootstrap-growing", {
    growing |>
      group_by(truth, growth, method) |>
      group_modify(function(rows, keys) slope_of_medians(rows)) |>
      ungroup() |>
      mutate(
        item = ifelse(truth == "null", "3b.8", "3b.9"),
        quantity = "slope of median log BF against log n",
        expected = "positive",
        delta = direction_delta(estimate, relative_delta_3b, direction_delta_3b)
      ) |>
      with_targets()
  })
}

growing_medians_3b <- function(growing) {
  growing |>
    group_by(truth, growth, n, method) |>
    summarise(median_log_bf = median(log_bf), .groups = "drop")
}


## ---- The results file ---------------------------------------------------------

## `shortfall`: by how much bridge sampling's reported error falls short of
## the real spread (study 3a, part B); 3b.1 and 3b.2 are about one given
## dataset each, so their errors are the reported ones times that factor
save_results_3b <- function(formula_checks, rates, growing, intrinsic_draws, shortfall) {
  formula_checks$bridge_error <- formula_checks$error
  formula_checks$error <- formula_checks$error * shortfall

  # 3b.1 predictive matching: expected 0
  at_minimum <- subset(formula_checks, item == "3b.1")
  matching <- expectation_rows(
    prior = at_minimum$prior,
    k = at_minimum$k,
    n = at_minimum$n,
    divergent = at_minimum$divergent,
    bridge_error = at_minimum$bridge_error,
    quantity = "3b.1 log BF at n = k + 1",
    expected = 0,
    estimate = at_minimum$log_bf,
    mcse = at_minimum$error
  )

  # 3b.2 information consistency: expected to increase
  information <- subset(formula_checks, item == "3b.2")
  information <- information[
    , c("prior", "r_squared", "log_bf", "error", "bridge_error", "divergent")
  ]
  information$expected <- "increases with R-squared"

  # 3b.3-3b.5 the rates
  slopes <- rate_slopes_3b(rates)

  # Stan's divergent transitions behind the three priors' Bayes factors
  divergences <- bind_rows(
    select(formula_checks, prior, divergent),
    select(rates, prior, divergent)
  ) |>
    group_by(prior) |>
    summarise(
      bayes_factors = n(),
      with_divergences = sum(divergent > 0),
      divergent_transitions = sum(divergent)
    )

  save_study_results_file("3b", list(
    predictive_matching = matching,
    information_consistency = information,
    rates = slopes,
    beyond_reference = rbind(beyond_reference(matching), beyond_reference(slopes)),
    alternative = alternative_3b(rates),
    alternative_slopes = alternative_slopes_3b(rates),
    intrinsic_consistency = intrinsic_3b(intrinsic_draws),
    growing_slopes = growing_slopes_3b(growing),
    growing_medians = growing_medians_3b(growing),
    divergences = divergences
  ))
}
