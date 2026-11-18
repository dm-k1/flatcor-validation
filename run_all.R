## =====================================================================
##  Run the whole validation: the pipeline (_targets.R), then the figures
##
##    Rscript run_all.R                  every study, then every figure
##    Rscript run_all.R 3a               one study (its targets end in _3a), then the figures
##    VALIDATION_SMOKE=1 Rscript run_all.R   the smoke test (results in cache/)
##    VALIDATION_FRESH=1 Rscript run_all.R   the whole validation from scratch
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
##  Each study writes results/studies/study_<cell>.rds; each
##  script in R/3_figures/ reads those files and writes its figures. The
##  notebook (validation.ipynb) presents the same results.
## =====================================================================
smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
fresh <- identical(Sys.getenv("VALIDATION_FRESH"), "1")
requested <- commandArgs(trailingOnly = TRUE)
store <- if (smoke) file.path("cache", "targets_smoke") else "_targets"
results_root <- if (smoke) file.path("cache", "smoke_results") else "results"
figure_folder <- if (smoke) file.path("cache", "smoke_figures") else file.path("results", "figures")

## ---- From scratch: always for the smoke test, on request for a real run -------------
if (fresh && length(requested) > 0) {
  stop(
    "VALIDATION_FRESH=1 starts the whole validation from scratch: run it without ",
    "naming studies (to redo one study, use run_studies(fresh = TRUE) in the notebook)."
  )
}
if (smoke || fresh) {
  unlink(c(store, file.path(results_root, "studies"), figure_folder), recursive = TRUE)
}

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
  design_scripts
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
    study = "1b", tables = "pooled", by = "prior",
    cells = priors_1b, first_round = first_round_1b, batch_size = batch_size_1b
  ),
  "2b" = list(
    study = "2b", tables = "estimates", by = c("k", "n", "rho", "others"),
    cells = cells_2b, first_round = first_round_2b, batch_size = batch_size_2b
  ),
  "3a" = list(
    study = "3a", tables = "calibration_line", by = c("prior", "k", "n"),
    cells = cells_3a, first_round = first_round_3a, batch_size = batch_size_3a
  ),
  "3a-repeats" = list(
    study = "3a", tables = "repeated_runs", by = c("prior", "k", "n"),
    cells = cells_3a, first_round = first_round_repeats_3a, batch_size = batch_size_repeats_3a
  ),
  "3b-rates" = list(
    study = "3b", tables = c("rates", "alternative_slopes"),
    by = c("item", "setting", "q", "prior"),
    cells = rate_cells_3b, first_round = first_round_rates_3b, batch_size = batch_size_3b
  ),
  "3b-growing" = list(
    study = "3b", tables = "growing_slopes", by = c("truth", "growth", "method"),
    cells = growing_cells_3b, first_round = first_round_growing_3b, batch_size = batch_size_3b
  )
)
if (length(requested) > 0) {
  in_rounds <- Filter(function(design) design$study %in% requested, in_rounds)
}
budget_core_hours <- 2500 # about 20 euros (docs/working/run_time_estimate.md)
last_round <- if (smoke) 2 else Inf

## the core hours so far: every piece's seconds, times the 4 chains it runs
core_hours_so_far <- function() {
  pieces <- targets::tar_meta(fields = c("type", "seconds"), store = store)
  sum(pieces$seconds[pieces$type %in% c("stem", "branch")], na.rm = TRUE) * 4 / 3600
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
  extended <- character(0)
  for (name in names(in_rounds)) {
    design <- in_rounds[[name]]
    plan <- read_plan(name, store, design$cells, design$first_round)
    file <- file.path(results_root, "studies", sprintf("study_%s.rds", design$study))
    tables <- readRDS(file)$tables[design$tables]
    columns <- c(design$by, "mcse", "target_mcse", "reached")
    table <- do.call(rbind, lapply(tables, function(t) as.data.frame(t)[columns]))
    plan_next <- next_plan(plan, table, design$by, design$batch_size, one_batch = smoke)
    if (!identical(plan_next$datasets, plan$datasets)) {
      write_plan(plan_next, name, store)
      extended <- c(extended, sprintf(
        "%s (%d of %d cells)", name, sum(plan_next$datasets > plan$datasets), nrow(plan)
      ))
    }
  }

  core_hours <- core_hours_so_far()
  message(sprintf(
    "Round %d done, %.1f core hours so far. %s", round, core_hours,
    if (length(extended) == 0) "Every cell has reached its targets." else
      paste("Next round for:", paste(extended, collapse = ", "))
  ))
  if (length(extended) == 0 || round >= last_round) break
  if (core_hours >= budget_core_hours) {
    message("The budget of ", budget_core_hours, " core hours is reached: stopping here. ",
            "Cells short of a target are marked in their study's results (column reached).")
    break
  }
  round <- round + 1
}

## ---- The figures, each script in its own R process ----------------------------------
results_folder <- file.path(results_root, "studies")
failed_figures <- character(0)
for (script in list.files(file.path("R", "3_figures"), pattern = "[.]R$", full.names = TRUE)) {
  study <- sub("^([0-9][a-z]).*", "\\1", basename(script))
  if (!file.exists(file.path(results_folder, sprintf("study_%s.rds", study)))) next
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  message(sprintf("%s: %s", script, if (status == 0) "done" else "stopped (see above)"))
  if (status != 0) failed_figures <- c(failed_figures, script)
}

## ---- Which results there are -------------------------------------------------------
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
check_results(results_folder, studies_run, flatcor_commit, failed_figures, smoke)

cat(sprintf("Run took %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
