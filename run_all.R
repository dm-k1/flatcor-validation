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

## ---- The studies -------------------------------------------------------------------
started <- Sys.time()
if (length(requested) == 0) {
  targets::tar_make(store = store)
} else {
  targets::tar_make(names = tidyselect::ends_with(!!paste0("_", requested)), store = store)
}

## ---- The figures, each script in its own R process ----------------------------------
results_folder <- file.path(results_root, "studies")
for (script in list.files(file.path("R", "3_figures"), pattern = "[.]R$", full.names = TRUE)) {
  study <- sub("^([0-9][a-z]).*", "\\1", basename(script))
  if (!file.exists(file.path(results_folder, sprintf("study_%s.rds", study)))) next
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  message(sprintf("%s: %s", script, if (status == 0) "done" else "stopped (see above)"))
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
cat(sprintf("Run took %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
