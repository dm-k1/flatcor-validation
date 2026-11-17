## =====================================================================
##  Helpers of the scripts behind the docs (R/docs/motivating_simulation.R):
##  parallel work with an independent random stream per item, saving figures,
##  and the simulation's own pass/fail checks (M0 to M4). Not part of the
##  validation, which runs as a pipeline (_targets.R).
## =====================================================================

## ---- Text encoding ------------------------------------------------------------------
## A cluster job or a Jupyter kernel may start R in an ASCII locale; the
## results contain characters such as ² and names such as García-Donato, so
## use UTF-8 (the results files then read
## back without warnings anywhere)
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (locale in c("C.UTF-8", "en_US.UTF-8", "en_GB.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break
  }
}

## ---- Parallel ----------------------------------------------------------------
## parallel_map(items, work) runs work(item) for every item on all cores but one.
## It forks the R session (fast, shares everything) where that is allowed: under
## Rscript on macOS and Linux. Where forking is not allowed (Windows; Positron
## and RStudio sessions), it starts fresh R processes instead and copies the
## global objects and attached packages into them.
## VALIDATION_CORES overrides it (a cloud machine or a cluster job that should
## use exactly the cores it was given).
number_of_cores <- as.integer(Sys.getenv(
  "VALIDATION_CORES",
  max(1, parallel::detectCores() - 1)
))
options(mc.cores = number_of_cores)

forking_allowed <- .Platform$OS.type != "windows" &&
  !inherits(try(parallel::mclapply(1:2, identity, mc.cores = 2), silent = TRUE), "try-error")

## Each item gets its own independent random stream (L'Ecuyer-CMRG, as in
## parallel::nextRNGStream), started from the caller's stream: a set.seed()
## before the call fixes them all, whatever the number of cores, so the same
## script gives the same numbers on any machine. (Seeding each item with
## set.seed(base + i) instead is not safe: consecutive seeds can give draws
## that are not quite uniform across items.)
parallel_map <- function(items, work, cores = number_of_cores) {
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
  ## an error in one item stops the run with that error, instead of coming
  ## back as a value that fails later
  stop_on_failure <- function(results) {
    failed <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      stop("parallel_map(): item ", which(failed)[1], " failed: ",
        conditionMessage(attr(results[[which(failed)[1]]], "condition")),
        call. = FALSE
      )
    }
    results
  }
  if (forking_allowed) {
    return(stop_on_failure(parallel::mclapply(seq_along(items), seeded_work, mc.cores = cores)))
  }
  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, setwd, getwd())
  ## attach the same packages; the global objects travel along below
  attached_packages <- .packages()
  parallel::clusterCall(cluster, function(packages) {
    for (package in rev(packages)) library(package, character.only = TRUE)
    NULL
  }, attached_packages)
  parallel::clusterExport(cluster, ls(globalenv()), envir = globalenv())
  parallel::parLapply(cluster, seq_along(items), seeded_work)
}

## ---- Figures -------------------------------------------------------------------
## save_figure() writes to figures_folder (the calling script sets it)
figures_folder <- file.path("docs", "figures")
## A white background with a light grid
ggplot2::theme_set(ggplot2::theme_bw())
## Saves a ggplot to <figures_folder>/<name>.png, and also shows it when
## running interactively. (Under Rscript a printed plot would otherwise end up
## in Rplots.pdf in the project folder.)
save_figure <- function(figure, name, width = 10, height = 6) {
  file <- file.path(figures_folder, paste0(name, ".png"))
  ggplot2::ggsave(file, figure, width = width, height = height, dpi = 150)
  if (interactive()) print(figure)
  invisible(figure)
}

## ---- The simulation's own checks ---------------------------------------------------
## record_check() records one check and prints PASS or FAIL; report_checks()
## prints them and stops with an error if one failed.
check_results <- data.frame(check = character(), passed = logical(), details = character())
record_check <- function(name, passed, details = "") {
  passed <- isTRUE(passed)
  # nolint next: assignment_linter. (the checks are collected in the script's environment)
  check_results[nrow(check_results) + 1, ] <<- list(name, passed, details)
  cat(sprintf("[%s] %-55s %s\n", if (passed) "PASS" else "FAIL", name, details))
}
report_checks <- function() {
  cat("\n")
  print(check_results, right = FALSE, row.names = FALSE)
  number_failed <- sum(!check_results$passed)
  if (number_failed > 0) {
    stop(sprintf("%d of %d checks FAILED.", number_failed, nrow(check_results)))
  }
  cat(sprintf("\nAll %d checks passed.\n", nrow(check_results)))
}
