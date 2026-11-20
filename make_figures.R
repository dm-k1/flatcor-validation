## =====================================================================
##  The figures, from the results files (after run_all.R)
##
##    Rscript make_figures.R                     every study that has a results file
##    Rscript make_figures.R 3a                  one study
##    VALIDATION_SMOKE=1 Rscript make_figures.R  the smoke test's figures (cache/)
##    VALIDATION_DRY=1 Rscript make_figures.R    the dry run's figures (cache/)
##
##  run_all.R only simulates, so a cloud machine spends its time on the
##  studies alone. The figures are drawn here, on any machine: after a run,
##  or after downloading its results (results/studies/study_<cell>.rds).
##  Each script in R/3_figures/ reads one study's results file, computes
##  nothing, and writes its figures to results/figures/ (cache/ for the
##  smoke test and the dry run), each in its own R process. The notebook
##  (validation.ipynb) runs this after the studies, and so does the smoke
##  test in CI.
##
##  A script that stops is an error here, after the others have run.
## =====================================================================
dry <- identical(Sys.getenv("VALIDATION_DRY"), "1")
# the dry run has the smoke test's design (set here, so the figure scripts inherit it)
if (dry) Sys.setenv(VALIDATION_SMOKE = "1")
smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
requested <- commandArgs(trailingOnly = TRUE)

## the same folders as run_all.R and R/shared/figures.R
results_root <- if (dry) {
  file.path("cache", "dry_results")
} else if (smoke) {
  file.path("cache", "smoke_results")
} else {
  "results"
}
results_folder <- file.path(results_root, "studies")

scripts <- list.files(file.path("R", "3_figures"), pattern = "[.]R$", full.names = TRUE)
failed <- character(0)
drawn <- 0
for (script in scripts) {
  study <- sub("^([0-9][a-z]).*", "\\1", basename(script))
  if (length(requested) > 0 && !(study %in% requested)) next
  if (!file.exists(file.path(results_folder, sprintf("study_%s.rds", study)))) {
    message(sprintf("%s: no results file for %s in %s, skipped", script, study, results_folder))
    next
  }
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  message(sprintf("%s: %s", script, if (status == 0) "done" else "stopped (see above)"))
  if (status == 0) drawn <- drawn + 1 else failed <- c(failed, script)
}

if (length(failed) > 0) {
  stop("Figure script(s) stopped: ", paste(failed, collapse = ", "), call. = FALSE)
}
message(sprintf("Figures: %d script(s) done.", drawn))
