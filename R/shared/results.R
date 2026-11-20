## =====================================================================
##  The results file of a study: results/studies/study_<cell>.rds
##  (cache/smoke_results/studies/ for the smoke test, cache/timing/results/
##  for a timing run, so neither mixes with real results), with the study's
##  tables and the setup of the run. The notebook and the figure scripts
##  read only these files.
## =====================================================================

## Writes the file and returns its path (a targets "file" target)
save_study_results_file <- function(study, tables) {
  smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
  timing <- !is.na(suppressWarnings(as.numeric(Sys.getenv("VALIDATION_TIMING", NA))))
  root <- if (smoke) {
    file.path("cache", "smoke_results")
  } else if (timing) {
    file.path("cache", "timing", "results")
  } else {
    "results"
  }
  folder <- file.path(root, "studies")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(folder, sprintf("study_%s.rds", study))

  # where, when and with what the study ran
  or_missing <- function(expression) {
    tryCatch(expression, error = function(e) NA_character_, warning = function(w) NA_character_)
  }
  commit <- or_missing(system2(
    "git", c("rev-parse", "--short", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ))
  setup <- list(
    finished = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    commit = if (length(commit) == 1) commit else NA_character_,
    R = R.version.string,
    flatcor = as.character(utils::packageVersion("flatcor")),
    flatcor_commit = utils::packageDescription("flatcor")$RemoteSha,
    CmdStan = or_missing(as.character(cmdstanr::cmdstan_version())),
    machine = paste(Sys.info()[["sysname"]], Sys.info()[["machine"]]),
    root_seed = root_seed, # every random number follows from it (R/shared/random_streams.R)
    smoke = smoke
  )

  saveRDS(list(study = study, tables = tables, setup = setup), file)
  file
}
