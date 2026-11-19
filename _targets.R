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

## flatcor compiles its Stan model into R's user cache folder: keep that in
## cache/, writable everywhere (set before the workers start, so they inherit it)
if (!nzchar(Sys.getenv("R_USER_CACHE_DIR"))) {
  Sys.setenv(R_USER_CACHE_DIR = file.path(getwd(), "cache", "r_user_cache"))
}

## Linux: every Stan model loads CmdStan's TBB library (libtbb.so.2) when it
## starts, and compiling a model can make CmdStan rebuild that library in
## place. With several workers, one worker's compile then removes it under a
## model another worker is starting ("libtbb.so.2: cannot open shared object
## file", seen with 7 workers on AWS). A copy made here, before any worker
## starts, and put first on the library path the workers inherit, is never
## touched by such a rebuild.
if (Sys.info()[["sysname"]] == "Linux") {
  tbb_copy <- file.path(getwd(), "cache", "tbb")
  dir.create(tbb_copy, recursive = TRUE, showWarnings = FALSE)
  tbb_libraries <- list.files(
    file.path(cmdstanr::cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb"),
    pattern = "^libtbb.*[.]so[.]2$", full.names = TRUE
  )
  file.copy(tbb_libraries, tbb_copy, overwrite = TRUE)
  library_path <- Sys.getenv("LD_LIBRARY_PATH")
  if (!startsWith(library_path, tbb_copy)) {
    library_path <- paste(c(tbb_copy, library_path[nzchar(library_path)]), collapse = ":")
    Sys.setenv(LD_LIBRARY_PATH = library_path)
  }
}

## One thread per process for linear algebra (BLAS, OpenMP). Each piece of
## work already keeps 4 cores busy with Stan's chains, and the workers are
## sized for that; a BLAS that starts a thread per core in every worker
## overloads the machine instead (a load of about 200 on 64 cores, seen on
## AWS). Set before the workers start, so they inherit it; a value set from
## outside (a cluster job, say) is kept.
local({
  for (variable in c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")) {
    if (!nzchar(Sys.getenv(variable))) do.call(Sys.setenv, stats::setNames(list("1"), variable))
  }
})

## Each piece of work runs Stan's 4 chains in parallel, so a quarter of the
## cores work on pieces at the same time, each in its own R process. The
## cores are those given to this run: on a cluster, the job's allocation
## (parallelly reads SLURM, PBS and cgroups), not all cores of the node.
## The dry run samples nothing, so it runs every step in this one process:
## starting workers would take longer than the steps.
tar_option_set(
  packages = c("flatcor", "dplyr"),
  controller = if (dry_run) NULL else crew_controller_local(
    workers = max(1, (parallelly::availableCores() - 1) %/% 4),
    # a log per worker (its R output, including Stan's own messages), for
    # finding out why a piece of work failed
    options_local = crew_options_local(log_directory = file.path("cache", "logs", "workers"))
  )
)

## The shared code, then each study's design and code
studies <- c(
  "1a_prior_sampling", "1b_invariance", "2a_posterior",
  "2b_estimation", "3a_bayes_factors", "3b_default_criteria"
)
tar_source(c(
  "R/shared/random_streams.R", "R/shared/pipeline.R", "R/shared/shared_stan_build.R",
  "R/shared/expectations.R",
  "R/shared/results.R", "R/shared/comparison_methods.R", "R/shared/rounds.R",
  list.files(file.path("R", "1_data"), pattern = "[.]R$", full.names = TRUE),
  file.path("R", "2_studies", studies, "design.R"),
  file.path("R", "2_studies", studies, "code.R"),
  # last, and only to record or replay the expensive calls: it replaces a few
  # of the functions above (R/dry_run/replay.R)
  if (dry_run || identical(Sys.getenv("VALIDATION_RECORD"), "1")) "R/dry_run/replay.R"
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
  # the plan and its pieces, as in 2b (R/shared/rounds.R)
  tar_target(
    plan_1b,
    read_plan("1b", targets::tar_path_store(), priors_1b, first_round_1b),
    cue = tar_cue(mode = "always")
  ),
  tar_target(pieces_1b, pieces_of_plan(plan_1b, batch_size_1b)),
  tar_target(
    changes_1b,
    invariance_batch_1b(pieces_1b["prior"], pieces_1b$from:pieces_1b$to),
    pattern = map(pieces_1b)
  ),
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
  tar_target(results_2a, save_results_2a(sbc_2a), format = "file"),


  ## ---- 2b. How close are the estimates to a known truth? ----
  # the plan: per cell, the datasets of the rounds so far (R/shared/rounds.R);
  # read every time, but what follows reruns only if it changed
  tar_target(
    plan_2b,
    read_plan("2b", targets::tar_path_store(), cells_2b, first_round_2b),
    cue = tar_cue(mode = "always")
  ),
  # one piece of work per cell and batch; a new round only adds pieces
  tar_target(pieces_2b, pieces_of_plan(plan_2b, batch_size_2b)),
  tar_target(
    fits_2b,
    estimates_batch_2b(pieces_2b[names(cells_2b)], pieces_2b$from:pieces_2b$to, stan_model),
    pattern = map(pieces_2b)
  ),
  tar_target(results_2b, save_results_2b(fits_2b), format = "file"),


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
    calibration_batch_3a(pieces_3a[names(cells_3a)], pieces_3a$from:pieces_3a$to),
    pattern = map(pieces_3a)
  ),
  # the data-averaged posterior decides; the bins and the line describe
  tar_target(averaged_table_3a, data_averaged_posterior_3a(calibration_runs_3a)),
  tar_target(bins_table_3a, calibration_bins_3a(calibration_runs_3a)),
  tar_target(line_table_3a, calibration_line_3a(calibration_runs_3a)),

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
    repeats_3a(
      pieces_repeats_3a[names(cells_3a)],
      pieces_repeats_3a$from:pieces_repeats_3a$to,
      repeat_runs_3a
    ),
    pattern = map(pieces_repeats_3a)
  ),
  tar_target(spread_table_3a, repeat_spread_3a(repeated_3a, repeat_runs_3a)),
  # the factor by which the reported error falls short, over every cell:
  # used by 1b and 3b for quantities about one given dataset
  tar_target(shortfall_3a, bridge_shortfall_3a(repeated_3a)),

  # also reported, and the results file
  tar_target(diagnostics_table_3a, diagnostics_3a(calibration_runs_3a)),
  tar_target(
    results_3a,
    save_results_3a(
      averaged_table_3a, bins_table_3a, line_table_3a, spread_table_3a, shortfall_3a,
      diagnostics_table_3a, calibration_runs_3a, repeated_3a
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
    batch_3b(
      pieces_rates_3b[names(rate_cells_3b)],
      pieces_rates_3b$from:pieces_rates_3b$to,
      rate_bayes_factor_3b
    ),
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
    batch_3b(
      pieces_growing_3b[names(growing_cells_3b)],
      pieces_growing_3b$from:pieces_growing_3b$to,
      growing_bayes_factor_3b
    ),
    pattern = map(pieces_growing_3b)
  ),

  # the results file (3b.7 is computed in it)
  tar_target(
    results_3b,
    save_results_3b(formula_3b, rates_3b, growing_3b, shortfall_3a),
    format = "file"
  )
)
