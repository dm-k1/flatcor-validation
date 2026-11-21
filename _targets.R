## =====================================================================
##  The validation as one pipeline (targets). This file is the "run all".
##
##    targets::tar_make()                          run everything that is not done yet
##    targets::tar_make(names = ends_with("_3a"))  one study
##    targets::tar_outdated()                      what tar_make() would (re)run
##    targets::tar_visnetwork()                    the steps and what depends on what
##
##  Each study writes its tables to results/studies/study_<cell>.rds; the
##  figures are made from those files afterwards by the scripts in
##  R/3_figures/ (Rscript run_all.R runs the studies' rounds through a queue,
##  cell by cell, then this pipeline once; then
##  Rscript make_figures.R the figures). Every finished piece is saved in _targets/
##  as soon as it is done, so a run that stops (a laptop closed, Colab or a
##  cloud machine shutting down) continues where it was.
##
##  Smoke test (a corner of each design, to test the code; its own store and
##  results folder, so it never mixes with a real run):
##    Sys.setenv(VALIDATION_SMOKE = "1"); targets::tar_make(store = "cache/targets_smoke")
##
##  Dry run (the smoke test's design with Stan and bridge sampling replayed
##  from a recorded smoke test, in one R process: tests the plumbing in a
##  minute or two; R/dry_run/replay.R):
##    Sys.setenv(VALIDATION_DRY = "1"); targets::tar_make(store = "cache/targets_dry")
##
##  Each study is a folder in R/2_studies/: design.R (the design) and code.R (the
##  steps). What each tests, expects and why: docs/validation_tests.md.
## =====================================================================
library(targets)
library(crew)


## ---- Setup ------------------------------------------------------------------------

## The dry run uses the smoke test's design (R/dry_run/replay.R)
dry_run <- identical(Sys.getenv("VALIDATION_DRY"), "1")
if (dry_run) Sys.setenv(VALIDATION_SMOKE = "1")

## Every worker's environment (the cache folder, a copy of the TBB library,
## one BLAS thread), set before any worker starts, so they inherit it, and
## the validation's code, in order (R/shared/workers.R)
source(file.path("R", "shared", "workers.R"))
tar_source(code_files())

## A worker per 4 cores (worker_count(), R/shared/workers.R). The dry run
## samples nothing, so it runs every step in this one process: starting
## workers would take longer than the steps.
tar_option_set(
  packages = c("flatcor", "dplyr"),
  # targets seeds R before each step from the step's name and this seed;
  # every draw that matters has its own seed (R/shared/random_streams.R),
  # but whatever is left follows the same root seed
  seed = root_seed,
  controller = if (dry_run) NULL else crew_controller_local(
    workers = worker_count(),
    # a log per worker (its R output, including Stan's own messages), for
    # finding out why a piece of work failed
    options_local = crew_options_local(log_directory = file.path("cache", "logs", "workers"))
  )
)


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
  # layer 1: the draws, the recipe's draws, each setting's run (R/shared/per_dataset.R)
  tar_target(
    per_dataset_draws_1a,
    write_per_dataset(draws_table_1a(samples_1a), "1a_prior_draws", gzip = TRUE),
    format = "file"
  ),
  tar_target(
    per_dataset_recipe_1a,
    write_per_dataset(
      recipe_table_1a(predictors_1a[["three_correlated"]]), "1a_recipe_draws", gzip = TRUE
    ),
    format = "file"
  ),
  tar_target(
    per_dataset_settings_1a,
    write_per_dataset(settings_table_1a(samples_1a), "1a_settings"),
    format = "file"
  ),
  tar_target(
    comparisons_1a,
    compare_1a(
      read_per_dataset(per_dataset_draws_1a), read_per_dataset(per_dataset_recipe_1a),
      read_per_dataset(per_dataset_settings_1a), predictors_1a
    )
  ),
  tar_target(results_1a, save_results_1a(comparisons_1a), format = "file"),


  ## ---- 1b. Is the Bayes factor invariant to units and coding? ----
  # the plan and its pieces, as in 2b (R/shared/rounds.R)
  tar_target(
    plan_1b,
    read_plan("1b", targets::tar_path_store(), priors_1b, first_round_1b),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_1b, pieces_of_plan(plan_1b, batch_size_1b)),
  tar_target(
    log_bfs_1b,
    piece_rows("1b", pieces_1b),
    pattern = map(pieces_1b)
  ),
  # layer 1: every log Bayes factor's row (R/shared/per_dataset.R)
  tar_target(
    per_dataset_1b,
    write_per_dataset(log_bfs_1b, "1b_log_bayes_factors"),
    format = "file"
  ),
  tar_target(changes_1b, changes_of_1b(read_per_dataset(per_dataset_1b))),
  # a change of one dataset is judged against bridge sampling's error as
  # measured in 3a, part B (shortfall_3a, below)
  tar_target(results_1b, save_results_1b(changes_1b, shortfall_3a), format = "file"),


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
  # layer 1: every fit's rows (R/shared/per_dataset.R)
  tar_target(
    per_dataset_fits_2a,
    write_per_dataset(bind_rows(lapply(sbc_2a, `[[`, "stats")), "2a_sbc_fits"),
    format = "file"
  ),
  tar_target(
    per_dataset_diagnostics_2a,
    write_per_dataset(bind_rows(lapply(sbc_2a, `[[`, "diagnostics")), "2a_sbc_diagnostics"),
    format = "file"
  ),
  tar_target(
    results_2a,
    save_results_2a(
      read_per_dataset(per_dataset_fits_2a), read_per_dataset(per_dataset_diagnostics_2a)
    ),
    format = "file"
  ),


  ## ---- 2b. How close are the estimates to a known truth? ----
  # the plan: per cell, the datasets of the rounds so far (R/shared/rounds.R);
  # read every time, but what follows reruns only if it changed
  tar_target(
    plan_2b,
    read_plan("2b", targets::tar_path_store(), cells_2b, first_round_2b),
    cue = tar_cue(mode = "always")
  ),
  # one piece of work per cell and batch; a new round only adds pieces. Each
  # piece's rows: the ones run_all.R's queue saved, or computed here
  # (piece_rows(), R/shared/rounds.R)
  tar_target(pieces_2b, pieces_of_plan(plan_2b, batch_size_2b)),
  tar_target(
    fits_2b,
    piece_rows("2b", pieces_2b, stan_model),
    pattern = map(pieces_2b)
  ),
  # layer 1: every fit's row (R/shared/per_dataset.R)
  tar_target(per_dataset_2b, write_per_dataset(fits_2b, "2b_fits"), format = "file"),
  tar_target(results_2b, save_results_2b(read_per_dataset(per_dataset_2b)), format = "file"),


  ## ---- 3a. Are the Bayes factors right? ----

  # part A, calibration: its plan and pieces, as in 2b (R/shared/rounds.R)
  tar_target(
    plan_3a,
    read_plan("3a", targets::tar_path_store(), cells_3a, first_round_3a),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_3a, pieces_of_plan(plan_3a, batch_size_3a)),
  tar_target(
    calibration_runs_3a,
    piece_rows("3a", pieces_3a),
    pattern = map(pieces_3a)
  ),
  # layer 1: every dataset's row (R/shared/per_dataset.R)
  tar_target(
    per_dataset_calibration_3a,
    write_per_dataset(calibration_runs_3a, "3a_calibration"),
    format = "file"
  ),
  tar_target(calibration_rows_3a, read_per_dataset(per_dataset_calibration_3a)),
  # the data-averaged posterior decides; the bins and the line describe
  tar_target(averaged_table_3a, data_averaged_posterior_3a(calibration_rows_3a)),
  tar_target(bins_table_3a, calibration_bins_3a(calibration_rows_3a)),
  tar_target(line_table_3a, calibration_line_3a(calibration_rows_3a)),

  # part B, bridge sampling's error: each dataset run repeat_runs_3a times;
  # its own plan and pieces
  tar_target(
    plan_repeats_3a,
    read_plan("3a-repeats", targets::tar_path_store(), cells_3a, first_round_repeats_3a),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_repeats_3a, pieces_of_plan(plan_repeats_3a, batch_size_repeats_3a)),
  tar_target(
    repeated_3a,
    piece_rows("3a-repeats", pieces_repeats_3a),
    pattern = map(pieces_repeats_3a)
  ),
  tar_target(
    per_dataset_repeats_3a,
    write_per_dataset(repeated_3a, "3a_repeated_runs"),
    format = "file"
  ),
  tar_target(repeated_rows_3a, read_per_dataset(per_dataset_repeats_3a)),
  tar_target(spread_table_3a, repeat_spread_3a(repeated_rows_3a, repeat_runs_3a)),
  # the factor by which the reported error falls short, over every cell:
  # used by 1b and 3b for quantities about one given dataset
  tar_target(shortfall_3a, bridge_shortfall_3a(repeated_rows_3a)),

  # also reported, and the results file
  tar_target(diagnostics_table_3a, diagnostics_3a(calibration_rows_3a)),
  tar_target(
    results_3a,
    save_results_3a(
      averaged_table_3a, bins_table_3a, line_table_3a, spread_table_3a, shortfall_3a,
      diagnostics_table_3a
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

  # 3b.3-3b.6: the rates; their plan and pieces, as in 2b (R/shared/rounds.R)
  tar_target(
    plan_rates_3b,
    read_plan("3b-rates", targets::tar_path_store(), rate_cells_3b, first_round_rates_3b),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_rates_3b, pieces_of_plan(plan_rates_3b, batch_size_3b)),
  tar_target(
    rates_3b,
    piece_rows("3b-rates", pieces_rates_3b),
    pattern = map(pieces_rates_3b)
  ),

  # 3b.8, 3b.9: k growing with n; their own plan and pieces
  tar_target(
    plan_growing_3b,
    read_plan("3b-growing", targets::tar_path_store(), growing_cells_3b, first_round_growing_3b),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_growing_3b, pieces_of_plan(plan_growing_3b, batch_size_3b)),
  tar_target(
    growing_3b,
    piece_rows("3b-growing", pieces_growing_3b),
    pattern = map(pieces_growing_3b)
  ),

  # 3b.7: the prior's draws with the sample and the population covariance
  tar_target(draws_intrinsic_3b, intrinsic_draws_of_3b()),

  # layer 1: every Bayes factor's row and the draws (R/shared/per_dataset.R)
  tar_target(per_dataset_formula_3b, write_per_dataset(formula_3b, "3b_formula"), format = "file"),
  tar_target(per_dataset_rates_3b, write_per_dataset(rates_3b, "3b_rates"), format = "file"),
  tar_target(
    per_dataset_growing_3b,
    write_per_dataset(growing_3b, "3b_growing_k", split_by = "growth"),
    format = "file"
  ),
  tar_target(
    per_dataset_intrinsic_3b,
    write_per_dataset(draws_intrinsic_3b, "3b_intrinsic_draws", gzip = TRUE),
    format = "file"
  ),

  # the results file, from the per-dataset files
  tar_target(
    results_3b,
    save_results_3b(
      read_per_dataset(per_dataset_formula_3b),
      read_per_dataset(per_dataset_rates_3b),
      read_per_dataset(per_dataset_growing_3b),
      read_per_dataset(per_dataset_intrinsic_3b),
      shortfall_3a
    ),
    format = "file"
  )
)
