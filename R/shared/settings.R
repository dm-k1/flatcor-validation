## =====================================================================
##  Shared settings for the validation scripts
##   - parallel: use all cores but one (leaves one free for the computer),
##               or VALIDATION_CORES if set
##   - caching:  simulation results are saved in cache/ (gitignored), one
##               folder per study and flatcor build, and reused
##   - figures:  saved to figures/ (tracked)
##   - checks:   record_check() and report_checks() for pass/fail summaries
## =====================================================================

## ---- Text encoding ------------------------------------------------------------------
## A cluster job or a Jupyter kernel may start R in an ASCII locale; the
## results contain characters such as ² and names such as García-Donato, so use UTF-8 (the results files then read
## back without warnings anywhere)
if (!isTRUE(l10n_info()$`UTF-8`))
  for (locale in c("C.UTF-8", "en_US.UTF-8", "en_GB.UTF-8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break

## ---- The package under test ------------------------------------------------------
## flatcor lives in its own repository (https://github.com/dm-k1/flatcor) and is
## installed at the version pinned in renv.lock. Some checks also call internal
## functions; they are taken from the package namespace here, in one place.
library(flatcor)
flatcor_internals <- c("sample_r_squared", "prior_shape", "log_sum_exp",
                       "log_likelihood_ratio_block")
for (internal in flatcor_internals)
  assign(internal, utils::getFromNamespace(internal, "flatcor"), envir = globalenv())

## ---- Parallel ----------------------------------------------------------------
## parallel_map(items, work) runs work(item) for every item on all cores but one.
## It forks the R session (fast, shares everything) where that is allowed: under
## Rscript on macOS and Linux. Where forking is not allowed (Windows; Positron
## and RStudio sessions), it starts fresh R processes instead and copies the
## global objects and attached packages into them.
## VALIDATION_CORES overrides it (a cloud machine or a cluster job that should
## use exactly the cores it was given).
number_of_cores <- as.integer(Sys.getenv("VALIDATION_CORES",
                                         max(1, parallel::detectCores() - 1)))
options(mc.cores = number_of_cores)

forking_allowed <- .Platform$OS.type != "windows" &&
  !inherits(try(parallel::mclapply(1:2, identity, mc.cores = 2), silent = TRUE), "try-error")

## Each item gets its own independent random stream (L'Ecuyer-CMRG, as in
## parallel::nextRNGStream), started from the caller's stream: a set.seed()
## before the call fixes them all, whatever the number of cores, so the same
## script gives the same numbers on any machine. (Seeding each item with
## set.seed(base + i) instead is not safe: consecutive seeds can give draws
## that are not quite uniform across items.)
parallel_map <- function(items, work) {
  start <- sample.int(.Machine$integer.max, 1)
  caller_kind <- RNGkind()[1]
  caller_seed <- get(".Random.seed", envir = globalenv())
  RNGkind("L'Ecuyer-CMRG")
  set.seed(start)
  streams <- vector("list", length(items))
  stream <- get(".Random.seed", envir = globalenv())
  for (index in seq_along(items)) {
    streams[[index]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  RNGkind(caller_kind)
  assign(".Random.seed", caller_seed, envir = globalenv())
  seeded_work <- function(index) {
    assign(".Random.seed", streams[[index]], envir = globalenv())
    work(items[[index]])
  }
  if (forking_allowed)
    return(parallel::mclapply(seq_along(items), seeded_work, mc.cores = number_of_cores))
  cluster <- parallel::makeCluster(number_of_cores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, setwd, getwd())
  ## attach the same packages (flatcor among them); the internal functions the
  ## checks use travel with the global objects below
  attached_packages <- .packages()
  parallel::clusterCall(cluster, function(packages) {
    for (package in rev(packages)) library(package, character.only = TRUE)
    NULL
  }, attached_packages)
  parallel::clusterExport(cluster, ls(globalenv()), envir = globalenv())
  parallel::parLapply(cluster, seq_along(items), seeded_work)
}

## ---- Folders -------------------------------------------------------------------
## Tracked: results/ (results/studies/ holds one study_NN.rds per study;
## run_all.R adds the run's summary and session information)
## and figures/. Everything that can be recomputed goes to cache/.
## Smoke mode (VALIDATION_SMOKE=1, see R/shared/precision.R) writes its
## results and figures to cache/ as well: they are tests of the code, not
## results.
## Timing mode (VALIDATION_TIMING=<fraction>, see cluster/time_study.R) runs
## the full design with full-length chains on that fraction of the datasets,
## to measure what the full run costs. Everything it writes, caches included,
## goes to cache/timing/<fraction>/, so it never mixes with a real run or with
## a timing run at another fraction.
smoke_mode      <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
timing_fraction <- suppressWarnings(as.numeric(Sys.getenv("VALIDATION_TIMING", NA)))
timing_mode     <- !smoke_mode && !is.na(timing_fraction)
cache_folder   <- if (timing_mode) file.path("cache", "timing", timing_fraction) else "cache"
results_folder <- if (smoke_mode) file.path(cache_folder, "smoke_results") else
                  if (timing_mode) file.path(cache_folder, "results") else "results"
figures_folder <- if (smoke_mode) file.path(cache_folder, "smoke_figures") else
                  if (timing_mode) file.path(cache_folder, "figures") else "figures"
studies_results_folder <- file.path(results_folder, "studies")
for (folder in c(studies_results_folder, figures_folder))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

## The build of flatcor under test: its commit when installed from GitHub (as
## renv.lock pins it), otherwise a hash of its functions and Stan model (a
## local install keeps the same version number while its code changes).
flatcor_build <- local({
  commit <- utils::packageDescription("flatcor")$RemoteSha
  if (!is.null(commit)) return(substr(commit, 1, 7))
  namespace <- asNamespace("flatcor")
  code <- lapply(sort(ls(namespace, all.names = TRUE)), function(name) {
    object <- get(name, envir = namespace)
    if (is.function(object)) deparse(object) else NULL
  })
  stan_model <- readLines(system.file("stan", "flatcor.stan", package = "flatcor"))
  paste0("local-", substr(digest::digest(list(code, stan_model)), 1, 7))
})

## A study's cache: cache/study_NN/flatcor_<build>/. The studies reuse a
## cached file whenever it exists, so each build of flatcor gets its own
## folder: results computed with another build are never reused. (A change to
## a study's design or number of datasets still needs its folder deleted.)
study_cache_folder <- function(study) {
  folder <- file.path(cache_folder, sprintf("study_%s%s", study, if (smoke_mode) "_smoke" else ""),
                      paste0("flatcor_", flatcor_build))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  folder
}

## ---- Figures -------------------------------------------------------------------
## A white background with a light grid, for every study's figures
ggplot2::theme_set(ggplot2::theme_bw())
## Saves a ggplot to <figures_folder>/<name>.png, and also shows it when
## running interactively. (Under Rscript a printed plot would otherwise end up
## in Rplots.pdf in the project folder.) The names are kept for the study's
## results file.
saved_figures <- character()
save_figure <- function(figure, name, width = 10, height = 6) {
  file <- file.path(figures_folder, paste0(name, ".png"))
  ggplot2::ggsave(file, figure, width = width, height = height, dpi = 150)
  saved_figures <<- union(saved_figures, file)
  if (interactive()) print(figure)
  invisible(figure)
}

## ---- Pass/fail checks ------------------------------------------------------------
## record_check() records one check and prints PASS or FAIL; report_checks()
## prints the summary and stops with an error if anything failed.
check_results <- data.frame(check = character(), passed = logical(), details = character())
record_check <- function(name, passed, details = "") {
  passed <- isTRUE(passed)
  check_results[nrow(check_results) + 1, ] <<- list(name, passed, details)
  cat(sprintf("[%s] %-55s %s\n", if (passed) "PASS" else "FAIL", name, details))
}
report_checks <- function() {
  cat("\n")
  print(check_results, right = FALSE, row.names = FALSE)
  number_failed <- sum(!check_results$passed)
  ## in smoke mode (R/shared/precision.R) the datasets are too few for the
  ## checks to mean anything: they are shown, but do not stop the run
  if (isTRUE(get0("smoke_mode", ifnotfound = FALSE)) || isTRUE(get0("timing_mode", ifnotfound = FALSE))) {
    cat(sprintf("\n%s mode: %d of %d checks would fail; with this few datasets that means nothing.\n",
                if (isTRUE(get0("smoke_mode", ifnotfound = FALSE))) "Smoke" else "Timing",
                number_failed, nrow(check_results)))
    return(invisible(check_results))
  }
  if (number_failed > 0) stop(sprintf("%d of %d checks FAILED.", number_failed, nrow(check_results)))
  cat(sprintf("\nAll %d checks passed.\n", nrow(check_results)))
}

## ---- The study's results file ----------------------------------------------------
## save_study_results("03", list(...)) writes results/studies/study_03.rds:
## the study's tables, its checks, the figures it saved and the setup of the
## run. The notebook (validation.ipynb) presents only these files.
## The studies may run on different machines (cluster/README.md), so each
## records its own: the processor, the cores it used, CmdStan and how long it took.
study_started <- Sys.time()
processor_name <- function() {
  name <- tryCatch(
    if (Sys.info()[["sysname"]] == "Darwin")
      system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE, stderr = FALSE) else
      sub(".*:\\s*", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
    error = function(e) NA_character_, warning = function(w) NA_character_)
  if (length(name) == 1) name else NA_character_
}
save_study_results <- function(study, tables) {
  git_commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                         error = function(e) NA_character_, warning = function(w) NA_character_)
  cmdstan <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_)
  saveRDS(list(study = study, tables = tables, checks = check_results, figures = saved_figures,
               setup = list(finished = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
                            commit = if (length(git_commit) == 1) git_commit else NA_character_,
                            R = R.version.string,
                            flatcor = as.character(utils::packageVersion("flatcor")),
                            flatcor_build = flatcor_build,
                            CmdStan = cmdstan,
                            machine = paste(Sys.info()[["sysname"]], Sys.info()[["machine"]]),
                            processor = processor_name(),
                            cores = number_of_cores,
                            minutes = round(as.numeric(difftime(Sys.time(), study_started, units = "mins")), 1),
                            smoke = smoke_mode,
                            timing = if (timing_mode) timing_fraction else NA)),
          file.path(studies_results_folder, sprintf("study_%s.rds", study)))
  session_folder <- file.path(results_folder, "session_info")
  dir.create(session_folder, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(session_folder, sprintf("study_%s.txt", study)))
}
