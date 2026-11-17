## =====================================================================
##  Study 2a: the code
##  What is tested and why: docs/validation_tests.md, study 2a.
##  The design: design.R. _targets.R puts these steps in order.
## =====================================================================


## ---- Simulation-based calibration of one batch -------------------------------------

## One SBC dataset: the true values of R-squared, b[1] and sigma, drawn from
## the prior (R/data/from_the_prior.R: the predictors, noise, intercept and
## sigma shared by the three priors), and the data for the Stan model
sbc_dataset_2a <- function(setting, replication) {
  data <- data_from_prior(
    setting$prior, setting$k, setting$n,
    x1 = "prior",
    purpose = "sbc",
    replication = replication,
    location_scale = "drawn"
  )
  truth <- list(
    R_squared = data$r_squared,
    `b[1]` = data$coefficients[1],
    sigma = data$noise_sd
  )
  stan_data <- list(
    N = setting$n,
    K = setting$k,
    X = data$predictors,
    y = data$outcome,
    prior_only = 0L,
    proper_stand_ins = 1L,
    r_squared_shape = prior_shape(setting$prior, setting$k)
  )

  true_values <- posterior::as_draws_rvars(lapply(truth, posterior::rvar))
  SBC::SBC_datasets(
    variables = posterior::as_draws_matrix(true_values),
    generated = list(stan_data)
  )
}

## Each fit: the posterior from flatcor's Stan model, and the rank of each
## true value among its posterior draws (the SBC package)
sbc_batch_2a <- function(setting, replications, model) {
  future::plan(future::sequential) # the pipeline runs the batches in parallel already

  datasets <- do.call(SBC::bind_datasets, lapply(replications, function(r) {
    sbc_dataset_2a(setting, r)
  }))

  backend <- SBC::SBC_backend_cmdstan_sample(
    model,
    chains = 1,
    iter_warmup = iterations_2a[["warmup"]],
    iter_sampling = iterations_2a[["sampling"]],
    seed = seed_of("2a", setting$prior, setting$k, replications[1]),
    refresh = 0,
    show_messages = FALSE
  )

  result <- SBC::compute_SBC(
    datasets, backend,
    keep_fits = FALSE,
    thin_ranks = iterations_2a[["sampling"]] / 100
  )

  list(
    stats = cbind(setting, as.data.frame(result$stats)),
    diagnostics = cbind(setting, as.data.frame(result$backend_diagnostics))
  )
}


## ---- The ranks against a uniform distribution -------------------------------------

## Scaled to 0-1: the mean (expected 1/2), the SD (expected that of a uniform
## rank) and a chi-square over 20 bins (expected 19, its SD under uniformity
## as the MCSE)
rank_summaries_2a <- function(stats) {
  summarise_ranks <- function(rows, keys) {
    scaled <- rows$rank / rows$max_rank[1]
    fits <- length(scaled)
    top <- rows$max_rank[1]
    uniform_sd <- sqrt(((top + 1)^2 - 1) / 12) / top

    bins <- table(cut(scaled, breaks = seq(0, 1, length.out = 21), include.lowest = TRUE))
    chi_square <- sum((bins - fits / 20)^2 / (fits / 20))

    tibble(
      quantity = c("mean scaled rank", "SD of the scaled rank", "chi-square over 20 rank bins"),
      expected = c(0.5, uniform_sd, 19),
      estimate = c(mean(scaled), sd(scaled), chi_square),
      mcse = c(
        uniform_sd / sqrt(fits),
        sqrt((1 / 80 - 1 / 144) / (4 / 12 * fits)),
        sqrt(2 * 19)
      )
    )
  }

  stats |>
    group_by(prior, k, n, variable) |>
    group_modify(summarise_ranks) |>
    ungroup() |>
    mutate(distance_in_mcse = (estimate - expected) / mcse)
}

## Stan's divergent transitions, per prior and setting
divergences_2a <- function(diagnostics) {
  diagnostics |>
    group_by(prior, k, n) |>
    summarise(
      fits = n(),
      fits_with_divergences = sum(n_divergent > 0),
      divergent_transitions = sum(n_divergent),
      .groups = "drop"
    )
}


## ---- The results file ---------------------------------------------------------

save_results_2a <- function(batches) {
  stats <- bind_rows(lapply(batches, `[[`, "stats"))
  diagnostics <- bind_rows(lapply(batches, `[[`, "diagnostics"))
  ranks <- rank_summaries_2a(stats)

  save_study_results_file("2a", list(
    ranks = ranks,
    beyond_reference = beyond_reference(ranks),
    divergences = divergences_2a(diagnostics),
    sbc_stats = stats
  ))
}
