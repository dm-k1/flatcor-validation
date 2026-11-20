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
##  until every cell reaches its targets or the run reaches its budget; the
##  smoke test runs two rounds.
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

## ---- The studies, in rounds --------------------------------------------------------
## The designs, for the studies that run in rounds. Each plan: its study's
## results file and the table(s) whose quantities must reach their targets,
## the columns that name a quantity's unit (a cell, or the cells of a slope
## over n), and its cells, first round and batch size
design_scripts <- list.files(
  file.path("R", "2_studies"),
  pattern = "^design[.]R$", recursive = TRUE, full.names = TRUE
)
for (script in c(
  file.path("R", "shared", c("random_streams.R", "pipeline.R", "rounds.R")),
  list.files(file.path("R", "1_data"), pattern = "[.]R$", full.names = TRUE),
  design_scripts,
  # last, and only to record or replay the expensive calls (R/dry_run/replay.R)
  if (dry || recording) file.path("R", "dry_run", "replay.R")
)) {
  source(script)
}

## ---- flatcor's models, compiled before any worker starts ---------------------------
## The first compile in a fresh CmdStan also builds the files every model
## shares (main.o, the TBB library). Done here, once, those files exist
## before the workers start, so a worker's own compile only writes its own
## files and never rebuilds a shared one under another worker. The lock in
## compile_flatcor_models() stays as a second guard. The cache folder is the
## one _targets.R gives the workers, so they find the model compiled here.
if (!nzchar(Sys.getenv("R_USER_CACHE_DIR"))) {
  Sys.setenv(R_USER_CACHE_DIR = file.path(getwd(), "cache", "r_user_cache"))
}
compile_flatcor_models()

in_rounds <- list(
  "1b" = list(
    study = "1b", tables = "pooled", by = "prior", pieces = "log_bfs_1b",
    cells = priors_1b, first_round = first_round_1b, batch_size = batch_size_1b
  ),
  "2b" = list(
    study = "2b", tables = "estimates", by = c("k", "n", "rho", "others"), pieces = "fits_2b",
    cells = cells_2b, first_round = first_round_2b, batch_size = batch_size_2b
  ),
  "3a" = list(
    study = "3a", tables = "data_averaged_posterior", by = c("prior", "k", "n"),
    pieces = "calibration_runs_3a",
    cells = cells_3a, first_round = first_round_3a, batch_size = batch_size_3a
  ),
  "3a-repeats" = list(
    study = "3a", tables = "repeated_runs", by = c("prior", "k", "n"), pieces = "repeated_3a",
    cells = cells_3a, first_round = first_round_repeats_3a, batch_size = batch_size_repeats_3a
  ),
  "3b-rates" = list(
    study = "3b", tables = c("rates", "alternative_slopes"),
    by = c("item", "setting", "q", "prior"), pieces = "rates_3b",
    cells = rate_cells_3b, first_round = first_round_rates_3b, batch_size = batch_size_3b
  ),
  "3b-growing" = list(
    study = "3b", tables = "growing_slopes", by = c("truth", "growth", "method"),
    pieces = "growing_3b",
    cells = growing_cells_3b, first_round = first_round_growing_3b, batch_size = batch_size_3b
  )
)
if (length(requested) > 0) {
  in_rounds <- Filter(function(design) design$study %in% requested, in_rounds)
}
## The budget: about 20 euros (docs/working/run_time_estimate.md). Raising it
## is a decision written in the log: VALIDATION_BUDGET_CORE_HOURS=3000
budget_core_hours <- as.numeric(Sys.getenv("VALIDATION_BUDGET_CORE_HOURS", "2500"))
message("Budget: ", budget_core_hours, " core hours")
last_round <- if (smoke) 2 else Inf

## the core hours so far: every piece's seconds, times the 4 chains it runs
core_hours_so_far <- function() {
  pieces <- targets::tar_meta(fields = c("type", "seconds"), store = store)
  sum(pieces$seconds[pieces$type %in% c("stem", "branch")], na.rm = TRUE) * 4 / 3600
}

## The core hours one dataset of a design has cost so far: the mean seconds
## of its pieces, times the 4 chains each runs, over the datasets in a piece
core_hours_per_dataset <- function(design) {
  meta <- targets::tar_meta(fields = c("type", "seconds"), store = store)
  pieces <- meta$seconds[meta$type == "branch" & startsWith(meta$name, paste0(design$pieces, "_"))]
  if (length(pieces) == 0) return(0)
  mean(pieces, na.rm = TRUE) * 4 / 3600 / design$batch_size
}

## The core hours the next round would cost, per design: the datasets it
## adds, times what one of its datasets has cost so far. `plans`,
## `plans_next`: per design, the plan so far and the plan of the next round.
round_cost <- function(plans, plans_next) {
  vapply(names(plans), function(name) {
    added <- sum(plans_next[[name]]$datasets - plans[[name]]$datasets)
    added * core_hours_per_dataset(in_rounds[[name]])
  }, numeric(1))
}

started <- Sys.time()
round <- 1
repeat {
  if (length(requested) == 0) {
    targets::tar_make(store = store)
  } else {
    targets::tar_make(names = tidyselect::ends_with(!!paste0("_", requested)), store = store)
  }

  # every cell short of a target gets the datasets it is estimated to lack
  plans <- list()
  plans_next <- list()
  for (name in names(in_rounds)) {
    design <- in_rounds[[name]]
    plans[[name]] <- read_plan(name, store, design$cells, design$first_round)
    file <- file.path(results_root, "studies", sprintf("study_%s.rds", design$study))
    tables <- readRDS(file)$tables[design$tables]
    columns <- c(design$by, "mcse", "target_mcse", "reached")
    table <- do.call(rbind, lapply(tables, function(t) as.data.frame(t)[columns]))
    plans_next[[name]] <- next_plan(
      plans[[name]], table, design$by, design$batch_size, one_batch = smoke
    )
  }
  short <- any(vapply(names(plans), function(name) {
    !identical(plans_next[[name]]$datasets, plans[[name]]$datasets)
  }, logical(1)))

  # a round that would take the run past its budget is not started: the run
  # stops before it and says what it would cost, and the researcher decides
  # (raise the budget and start again, or report the cells short of a
  # target as such); docs/validation_tests.md, "How many datasets"
  core_hours <- core_hours_so_far()
  cost <- round_cost(plans, plans_next)
  if (short && round < last_round && core_hours + sum(cost) > budget_core_hours) {
    message(sprintf(
      paste(
        "Round %d done, %.1f core hours so far. The next round would cost about %.1f core",
        "hours (%s), more than the %.1f left of the budget of %.1f: stopping before it.",
        "To run it, raise the budget (VALIDATION_BUDGET_CORE_HOURS) and start again; to",
        "stop here, the cells short of a target are marked in their study's results",
        "(column reached)."
      ),
      round, core_hours, sum(cost),
      paste(sprintf("%s %.1f", names(cost), cost)[cost > 0], collapse = ", "),
      max(0, budget_core_hours - core_hours), budget_core_hours
    ))
    break
  }

  extended <- character(0)
  for (name in names(plans)) {
    if (!identical(plans_next[[name]]$datasets, plans[[name]]$datasets)) {
      write_plan(plans_next[[name]], name, store)
      extended <- c(extended, sprintf(
        "%s (%d of %d cells)", name,
        sum(plans_next[[name]]$datasets > plans[[name]]$datasets), nrow(plans[[name]])
      ))
    }
  }

  status <- if (!short) {
    "Every cell has reached its targets."
  } else {
    paste("Next round for:", paste(extended, collapse = ", "))
  }
  message(sprintf("Round %d done, %.1f core hours so far. %s", round, core_hours, status))
  if (length(extended) == 0 || round >= last_round) break
  round <- round + 1
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
