## =====================================================================
##  Run the whole validation: the pipeline (_targets.R), in rounds
##
##    Rscript run_all.R                  every study
##    Rscript run_all.R 3a               one study (its targets end in _3a)
##    VALIDATION_SMOKE=1 Rscript run_all.R   the smoke test (results in cache/)
##    VALIDATION_SMOKE=1 VALIDATION_RECORD=1 Rscript run_all.R
##                                       the smoke test, recording the expensive calls
##    VALIDATION_DRY=1 Rscript run_all.R     the dry run: the smoke test's design with
##                                       those calls replayed, a minute or two (R/dry_run/replay.R)
##    VALIDATION_FRESH=1 Rscript run_all.R   the whole validation from scratch
##    VALIDATION_BUDGET_CORE_HOURS=3000 Rscript run_all.R   with a raised budget
##
##  The pipeline runs only what is not done yet: a run that stopped continues
##  where it was. It refuses to continue with another flatcor than the one
##  the run started with, so a results file never mixes two versions; start
##  from scratch then. VALIDATION_FRESH=1 deletes the store, the results and
##  the figures first. The smoke test always starts from scratch.
##  Studies that measure their numbers of datasets run in rounds
##  (R/shared/rounds.R; docs/validation_tests.md, "How many datasets"),
##  each cell on its own: a queue (R/shared/queue.R) runs their pieces of
##  work and gives a cell its next batch as soon as its own batch is done,
##  without waiting for the other cells, until every cell reaches its
##  targets or its next batch would take the run past the budget; in the
##  smoke test a cell gets at most one more batch. The parts with no rounds
##  (1a, 2a, 3b.1 and 3b.2) join the queue first and run alongside. Then the
##  pipeline runs once: every study's results, from the pieces. A study run
##  alone (Rscript run_all.R 1b) gives the same results as in a full run:
##  1b and 3b bring 3a part B along, whose measured bridge sampling error
##  they use.
##  Each study writes results/studies/study_<cell>.rds.
##
##  Only the simulation runs here, so a cloud machine does nothing else. The
##  figures come after, from those files: Rscript make_figures.R (with the
##  same VALIDATION_SMOKE or VALIDATION_DRY), on any machine, also after
##  downloading the results. The notebook (validation.ipynb) runs both, and
##  so does the smoke test in CI.
## =====================================================================
dry <- identical(Sys.getenv("VALIDATION_DRY"), "1")
recording <- identical(Sys.getenv("VALIDATION_RECORD"), "1")
if (dry && recording) stop("VALIDATION_DRY=1 replays recordings; it cannot also record them.")
if (recording && !identical(Sys.getenv("VALIDATION_SMOKE"), "1")) {
  stop("Recordings for the dry run come from a smoke test: VALIDATION_SMOKE=1 VALIDATION_RECORD=1")
}
# the dry run has the smoke test's design (set here, so the pipeline inherits it)
if (dry) Sys.setenv(VALIDATION_SMOKE = "1")
smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
fresh <- identical(Sys.getenv("VALIDATION_FRESH"), "1")
requested <- commandArgs(trailingOnly = TRUE)
store <- if (dry) {
  file.path("cache", "targets_dry")
} else if (smoke) {
  file.path("cache", "targets_smoke")
} else {
  "_targets"
}
results_root <- if (dry) {
  file.path("cache", "dry_results")
} else if (smoke) {
  file.path("cache", "smoke_results")
} else {
  "results"
}
# the figures of an earlier run, removed with its results when starting
# from scratch (make_figures.R draws new ones)
figure_folder <- if (dry) {
  file.path("cache", "dry_figures")
} else if (smoke) {
  file.path("cache", "smoke_figures")
} else {
  file.path("results", "figures")
}

## ---- From scratch: always for the smoke test, on request for a real run -------------
if (fresh && length(requested) > 0) {
  stop(
    "VALIDATION_FRESH=1 starts the whole validation from scratch: run it without ",
    "naming studies. To redo one study, mark its steps outdated first, ",
    "targets::tar_invalidate(names = tidyselect::ends_with(\"_3a\")), then Rscript run_all.R 3a."
  )
}
if (smoke || fresh) {
  unlink(c(store, file.path(results_root, "studies"), figure_folder), recursive = TRUE)
}
## A recording run replaces every recording, so none is left from an older
## design or flatcor (R/dry_run/replay.R)
if (recording) unlink(file.path("cache", "dry_recordings"), recursive = TRUE)

## ---- One flatcor per run ------------------------------------------------------------
## The store remembers the flatcor commit it was started with; continuing
## with another would mix two versions in one results file.
flatcor_commit <- utils::packageDescription("flatcor")$RemoteSha
if (is.null(flatcor_commit)) flatcor_commit <- "not installed from GitHub"
commit_file <- file.path(store, "flatcor_commit")
if (file.exists(commit_file)) {
  started_with <- readLines(commit_file)
  if (!identical(started_with, flatcor_commit)) {
    stop(
      "This run was started with flatcor ", started_with, "; installed now is ",
      flatcor_commit, ". Continuing would mix the two in the results. Install the ",
      "first again (renv::restore()), or start from scratch: VALIDATION_FRESH=1 Rscript run_all.R"
    )
  }
} else {
  dir.create(store, recursive = TRUE, showWarnings = FALSE)
  writeLines(flatcor_commit, commit_file)
}

## ---- The code ---------------------------------------------------------------------
## Every worker's environment and the validation's code (R/shared/workers.R),
## the same as the pipeline's (_targets.R), with the packages its steps
## load, and the queue (R/shared/queue.R)
suppressPackageStartupMessages({
  library(flatcor)
  library(dplyr)
})
source(file.path("R", "shared", "workers.R"))
for (script in c(code_files(), file.path("R", "shared", "queue.R"))) source(script)

## ---- flatcor's models, compiled before any worker starts ---------------------------
## The first compile in a fresh CmdStan also builds the files every model
## shares (main.o, the TBB library). Done here, once, those files exist
## before the workers start, so a worker's own compile only writes its own
## files and never rebuilds a shared one under another worker. The lock in
## compile_flatcor_models() stays as a second guard. The cache folder is the
## one R/shared/workers.R gives the workers, so they find the model compiled here.
compile_flatcor_models()

## ---- The studies that run in rounds: a queue, cell by cell (R/shared/queue.R) -------
## Per study (or part of one): its cells, first round and batch size; the
## columns that name a cell (one cell, or the cells of a slope over n); and
## `table`, the study's table(s) whose quantities must reach their targets,
## computed from a cell's per-dataset rows (the same functions its results
## file uses)
in_rounds <- list(
  "1b" = list(
    study = "1b", by = "prior",
    cells = priors_1b, first_round = first_round_1b, batch_size = batch_size_1b,
    table = function(rows) list(pooled_1b(changes_of_1b(rows)))
  ),
  "2b" = list(
    study = "2b", by = c("k", "n", "rho", "others"),
    cells = cells_2b, first_round = first_round_2b, batch_size = batch_size_2b,
    table = function(rows) list(results_table_2b(bias_and_coverage_2b(rows), prediction_2b(rows)))
  ),
  "3a" = list(
    study = "3a", by = c("prior", "k", "n"),
    cells = cells_3a, first_round = first_round_3a, batch_size = batch_size_3a,
    table = function(rows) list(data_averaged_posterior_3a(rows))
  ),
  "3a-repeats" = list(
    study = "3a", by = c("prior", "k", "n"),
    cells = cells_3a, first_round = first_round_repeats_3a, batch_size = batch_size_repeats_3a,
    table = function(rows) list(repeat_spread_3a(rows, repeat_runs_3a))
  ),
  "3b-rates" = list(
    study = "3b", by = c("item", "setting", "q", "prior"),
    cells = rate_cells_3b, first_round = first_round_rates_3b, batch_size = batch_size_3b,
    # a slope is either a rate (3b.3-3b.5) or the alternative's (3b.6)
    table = function(rows) {
      if (all(rows$setting == "alternative")) list(alternative_slopes_3b(rows))
      else list(rate_slopes_3b(rows))
    }
  ),
  "3b-growing" = list(
    study = "3b", by = c("truth", "growth", "method"),
    cells = growing_cells_3b, first_round = first_round_growing_3b, batch_size = batch_size_3b,
    table = function(rows) list(growing_slopes_3b(rows))
  )
)
for (name in names(in_rounds)) in_rounds[[name]]$name <- name

## The parts that run once, with no rounds: per part its study and its
## pieces (one setting of 1a; one batch of fits of a setting of 2a; one
## dataset of 3b.1 or 3b.2). They join the queue first, so they run
## alongside the rounds and nothing waits for them at the end.
run_once <- list(
  "1a" = list(study = "1a", pieces = settings_1a),
  "2a" = list(study = "2a", pieces = pieces_2a()),
  "3b-formula" = list(study = "3b", pieces = formula_checks_3b)
)

if (length(requested) > 0) {
  # 1b and 3b judge single datasets against bridge sampling's error as
  # measured by 3a part B (shortfall_3a in _targets.R), so part B runs to
  # its targets with them, as in a full run: a study run alone gives the
  # same results
  part_b <- any(c("1b", "3b") %in% requested)
  in_rounds <- Filter(function(design) {
    design$study %in% requested || (part_b && design$name == "3a-repeats")
  }, in_rounds)
  run_once <- Filter(function(part) part$study %in% requested, run_once)
}
## The budget: about 20 euros (docs/working/run_time_estimate.md). Raising it
## is a decision written in the log: VALIDATION_BUDGET_CORE_HOURS=3000
budget_core_hours <- as.numeric(Sys.getenv("VALIDATION_BUDGET_CORE_HOURS", "2500"))
message("Budget: ", budget_core_hours, " core hours")

## the core hours the pipeline's own steps have cost so far (a run that
## continues): every step's seconds, times the 4 chains it runs
pipeline_core_hours <- function() {
  if (!file.exists(file.path(store, "meta", "meta"))) return(0)
  steps <- targets::tar_meta(fields = c("type", "seconds"), store = store)
  sum(steps$seconds[steps$type %in% c("stem", "branch")], na.rm = TRUE) * 4 / 3600
}

started <- Sys.time()
run_rounds(
  in_rounds, lapply(run_once, `[[`, "pieces"), store, budget_core_hours,
  # the dry run samples nothing: every piece in this process
  workers = if (dry) 0 else worker_count(),
  spent_elsewhere = pipeline_core_hours(),
  one_extension = smoke
)

## ---- The pipeline: every study's results ---------------------------------------
## The pieces the queue ran are read from the store (piece_rows(),
## R/shared/rounds.R), not run again; what is left takes seconds
if (length(requested) == 0) {
  targets::tar_make(store = store)
} else {
  targets::tar_make(names = tidyselect::ends_with(!!paste0("_", requested)), store = store)
}

## ---- Which results there are -------------------------------------------------------
results_folder <- file.path(results_root, "studies")
for (file in sort(list.files(results_folder, pattern = "^study_.*[.]rds$", full.names = TRUE))) {
  setup <- readRDS(file)$setup
  cat(sprintf(
    "%s: finished %s, commit %s, flatcor %s (%s)%s\n", basename(file), setup$finished,
    setup$commit, setup$flatcor, substr(setup$flatcor_commit, 1, 7),
    if (isTRUE(setup$smoke)) " (smoke test)" else ""
  ))
}

## ---- Every output there, every number a number (R/shared/check_results.R) -----------
source(file.path("R", "shared", "check_results.R"))
studies_run <- if (length(requested) > 0) requested else
  sub("^([0-9][a-z]).*", "\\1", list.files(file.path("R", "2_studies")))
check_results(results_folder, studies_run, flatcor_commit, smoke)

cat(sprintf("Run took %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
