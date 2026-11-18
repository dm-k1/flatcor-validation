## =====================================================================
##  The validation as one pipeline (targets). This file is the "run all".
##
##    targets::tar_make()                          run everything that is not done yet
##    targets::tar_make(names = ends_with("_3a"))  one study
##    targets::tar_outdated()                      what tar_make() would (re)run
##    targets::tar_visnetwork()                    the steps and what depends on what
##
##  Each study writes its tables to results/studies/study_<cell>.rds; the
##  figures are made from those files by the scripts in R/3_figures/
##  (Rscript run_all.R does both). Every finished piece is saved in _targets/
##  as soon as it is done, so a run that stops (a laptop closed, Colab or a
##  cloud machine shutting down) continues where it was.
##
##  Smoke test (a corner of each design, to test the code; its own store and
##  results folder, so it never mixes with a real run):
##    Sys.setenv(VALIDATION_SMOKE = "1"); targets::tar_make(store = "cache/targets_smoke")
##
##  Each study is a folder in R/2_studies/: design.R (the design) and code.R (the
##  steps). What each tests, expects and why: docs/validation_tests.md.
## =====================================================================
library(targets)
library(crew)


## ---- Setup ------------------------------------------------------------------------

## flatcor compiles its Stan model into R's user cache folder: keep that in
## cache/, writable everywhere (set before the workers start, so they inherit it)
if (!nzchar(Sys.getenv("R_USER_CACHE_DIR"))) {
  Sys.setenv(R_USER_CACHE_DIR = file.path(getwd(), "cache", "r_user_cache"))
}

## Each piece of work runs Stan's 4 chains in parallel, so a quarter of the
## cores work on pieces at the same time, each in its own R process
tar_option_set(
  packages = c("flatcor", "dplyr"),
  controller = crew_controller_local(workers = max(1, (parallel::detectCores() - 1) %/% 4))
)

## The shared code, then each study's design and code
studies <- c(
  "1a_prior_sampling", "1b_invariance", "2a_posterior",
  "2b_estimation", "3a_bayes_factors", "3b_default_criteria"
)
tar_source(c(
  "R/shared/random_streams.R", "R/shared/pipeline.R", "R/shared/expectations.R",
  "R/shared/results.R", "R/shared/comparison_methods.R",
  list.files(file.path("R", "1_data"), pattern = "[.]R$", full.names = TRUE),
  file.path("R", "2_studies", studies, "design.R"),
  file.path("R", "2_studies", studies, "code.R")
))


## ---- The steps ------------------------------------------------------------------------

list(
  ## flatcor's Stan model, compiled once before the studies that sample it
  tar_target(stan_model, flatcor_stan_model()),


  ## ---- 1a. Does the software draw the prior it claims? ----
  tar_target(design_1a, settings_1a),
  tar_target(predictors_1a, data_prior_only_designs()),
  tar_target(
    samples_1a,
    sample_setting_1a(design_1a, predictors_1a[[design_1a$design]], stan_model),
    pattern = map(design_1a),
    iteration = "list"
  ),
  tar_target(
    comparisons_1a,
    compare_setting_1a(samples_1a, predictors_1a[[samples_1a$setting$design]]),
    pattern = map(samples_1a),
    iteration = "list"
  ),
  tar_target(results_1a, save_results_1a(comparisons_1a), format = "file"),


  ## ---- 1b. Is the Bayes factor invariant to units and coding? ----
  tar_target(design_1b, priors_1b),
  tar_target(
    batches_1b,
    split_into_batches(datasets_1b, batch_size_1b),
    iteration = "list"
  ),
  tar_target(
    changes_1b,
    invariance_batch_1b(design_1b, batches_1b),
    pattern = cross(design_1b, batches_1b)
  ),
  tar_target(results_1b, save_results_1b(changes_1b), format = "file"),


  ## ---- 2a. Is the posterior computed correctly? ----
  tar_target(design_2a, settings_2a),
  tar_target(
    batches_2a,
    split_into_batches(fits_2a, batch_size_2a),
    iteration = "list"
  ),
  tar_target(
    sbc_2a,
    sbc_batch_2a(design_2a, batches_2a, stan_model),
    pattern = cross(design_2a, batches_2a),
    iteration = "list"
  ),
  tar_target(results_2a, save_results_2a(sbc_2a), format = "file"),


  ## ---- 2b. How close are the estimates to a known truth? ----
  tar_target(design_2b, cells_2b),
  tar_target(
    batches_2b,
    split_into_batches(datasets_per_cell_2b, batch_size_2b),
    iteration = "list"
  ),
  tar_target(
    fits_2b,
    estimates_batch_2b(design_2b, batches_2b, stan_model),
    pattern = cross(design_2b, batches_2b)
  ),
  tar_target(results_2b, save_results_2b(fits_2b), format = "file"),


  ## ---- 3a. Are the Bayes factors right? ----
  tar_target(design_3a, cells_3a),

  # part A, calibration: every cell x every batch of datasets
  tar_target(
    batches_3a,
    split_into_batches(datasets_per_cell_3a, batch_size_3a),
    iteration = "list"
  ),
  tar_target(
    calibration_runs_3a,
    calibration_batch_3a(design_3a, batches_3a),
    pattern = cross(design_3a, batches_3a)
  ),
  tar_target(bins_table_3a, calibration_bins_3a(calibration_runs_3a)),
  tar_target(line_table_3a, calibration_line_3a(calibration_runs_3a)),

  # part B, bridge sampling's error: every cell, its datasets run repeat_runs_3a times
  tar_target(
    repeated_3a,
    repeats_3a(design_3a, seq_len(repeat_datasets_3a), repeat_runs_3a),
    pattern = map(design_3a)
  ),
  tar_target(spread_table_3a, repeat_spread_3a(repeated_3a, repeat_runs_3a)),

  # also reported, and the results file
  tar_target(diagnostics_table_3a, diagnostics_3a(calibration_runs_3a)),
  tar_target(
    results_3a,
    save_results_3a(
      bins_table_3a, line_table_3a, spread_table_3a, diagnostics_table_3a,
      calibration_runs_3a, repeated_3a
    ),
    format = "file"
  ),


  ## ---- 3b. Do the Bayes factors meet the default-prior criteria? ----

  # 3b.1, 3b.2: one dataset per row
  tar_target(design_formula_3b, formula_checks_3b),
  tar_target(
    formula_3b,
    formula_check_3b(design_formula_3b),
    pattern = map(design_formula_3b)
  ),

  # 3b.3-3b.6: the rates
  tar_target(design_rates_3b, rate_cells_3b),
  tar_target(
    rate_batches_3b,
    split_into_batches(rate_datasets_3b, batch_size_3b),
    iteration = "list"
  ),
  tar_target(
    rates_3b,
    batch_3b(design_rates_3b, rate_batches_3b, rate_bayes_factor_3b),
    pattern = cross(design_rates_3b, rate_batches_3b)
  ),

  # 3b.8, 3b.9: k growing with n
  tar_target(design_growing_3b, growing_cells_3b),
  tar_target(
    growing_batches_3b,
    split_into_batches(growing_datasets_3b, batch_size_3b),
    iteration = "list"
  ),
  tar_target(
    growing_3b,
    batch_3b(design_growing_3b, growing_batches_3b, growing_bayes_factor_3b),
    pattern = cross(design_growing_3b, growing_batches_3b)
  ),

  # the results file (3b.7 is computed in it)
  tar_target(results_3b, save_results_3b(formula_3b, rates_3b, growing_3b), format = "file")
)
