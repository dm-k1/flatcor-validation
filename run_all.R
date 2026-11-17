## =====================================================================
##  Run the whole validation: the pipeline (_targets.R), then the figures
##
##    Rscript run_all.R                  every study, then every figure
##    Rscript run_all.R 3a               one study (its targets end in _3a), then the figures
##    VALIDATION_SMOKE=1 Rscript run_all.R   the smoke test (results in cache/)
##
##  The pipeline runs only what is not done yet: a run that stopped continues
##  where it was. Each study writes results/studies/study_<cell>.rds; each
##  script in figures/ reads those files and writes its figures. The
##  notebook (validation.ipynb) presents the same results.
## =====================================================================
smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
requested <- commandArgs(trailingOnly = TRUE)
store <- if (smoke) file.path("cache", "targets_smoke") else "_targets"

## ---- The studies -------------------------------------------------------------------
started <- Sys.time()
if (length(requested) == 0) {
  targets::tar_make(store = store)
} else {
  targets::tar_make(names = tidyselect::ends_with(!!paste0("_", requested)), store = store)
}

## ---- The figures, each script in its own R process ----------------------------------
results_folder <- if (smoke) file.path("cache", "smoke_results") else "results"
results_folder <- file.path(results_folder, "studies")
for (script in list.files("figures", pattern = "[.]R$", full.names = TRUE)) {
  study <- sub("^([0-9][a-z]).*", "\\1", basename(script))
  if (!file.exists(file.path(results_folder, sprintf("study_%s.rds", study)))) next
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  message(sprintf("%s: %s", script, if (status == 0) "done" else "stopped (see above)"))
}

## ---- Which results there are -------------------------------------------------------
for (file in sort(list.files(results_folder, pattern = "^study_.*[.]rds$", full.names = TRUE))) {
  setup <- readRDS(file)$setup
  cat(sprintf(
    "%s: finished %s, commit %s, flatcor %s%s\n", basename(file), setup$finished, setup$commit,
    setup$flatcor, if (isTRUE(setup$smoke)) " (smoke test)" else ""
  ))
}
cat(sprintf("Run took %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
