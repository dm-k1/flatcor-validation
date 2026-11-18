## =====================================================================
##  The check at the end of every run (run_all.R): are all outputs there,
##  and is every number in them a number?
##
##    - a results file for every study that was run, all with the flatcor
##      the run used
##    - every table of every file: no NaN, no Inf, no NA except the blanks
##      listed below (each with its reason), and no empty table except
##      those listed below
##    - every figure script ran
##
##  The smoke test (and so every push to main, check.yml, and every image
##  build) stops with an error when something is wrong; a real run reports
##  it in its log, as its results are written and should not be lost.
## =====================================================================

## Blanks by design: the rows (when) where a column is NA on purpose
blanks_by_design <- list(
  list(
    study = "2a", table = "ranks", columns = c("delta", "target_mcse"),
    when = function(t) grepl("chi-square", t$quantity),
    reason = "the chi-square is an overall check, without a target of its own"
  ),
  list(
    study = "2b", table = "per_dataset", columns = "divergent",
    when = function(t) t$method %in% c("OLS", "JZS default", "Zellner-Siow"),
    reason = "these methods do not use Stan, so have no divergences to count"
  )
)

## Tables empty by design in the smoke test, whose corner of the design
## leaves their cells out
empty_in_smoke <- list(
  "3b" = c("alternative", "alternative_slopes") # no alternative setting in the smoke test
)

## The problems in one study's results file (character(0) for none)
check_results_file <- function(file, study, flatcor_commit, smoke) {
  if (!file.exists(file)) return(sprintf("%s: no results file", study))
  results <- readRDS(file)
  problems <- character(0)

  if (!identical(results$setup$flatcor_commit, flatcor_commit)) {
    problems <- c(problems, sprintf(
      "%s: made with flatcor %s, not %s",
      study, results$setup$flatcor_commit, flatcor_commit
    ))
  }

  for (name in names(results$tables)) {
    table <- as.data.frame(results$tables[[name]])
    where <- sprintf("%s, table %s", study, name)

    if (nrow(table) == 0) {
      if (!(smoke && name %in% empty_in_smoke[[study]])) {
        problems <- c(problems, sprintf("%s: empty", where))
      }
      next
    }

    # the blanks this table may have
    allowed <- matrix(FALSE, nrow(table), ncol(table), dimnames = list(NULL, names(table)))
    for (blank in blanks_by_design) {
      if (blank$study == study && blank$table == name) {
        allowed[blank$when(table), blank$columns] <- TRUE
      }
    }

    for (column in names(table)) {
      values <- table[[column]]
      is_nan <- if (is.numeric(values)) is.nan(values) else rep(FALSE, length(values))
      nan <- sum(is_nan)
      inf <- if (is.numeric(values)) sum(is.infinite(values)) else 0
      na <- sum(is.na(values) & !is_nan & !allowed[, column])
      if (nan + inf + na > 0) {
        problems <- c(problems, sprintf(
          "%s, column %s: %d NA, %d NaN, %d Inf", where, column, na, nan, inf
        ))
      }
    }
  }
  problems
}

## Checks the results of the studies that were run and the figure scripts;
## stops in the smoke test when anything is wrong
check_results <- function(results_folder, studies, flatcor_commit, failed_figures, smoke) {
  problems <- character(0)
  for (study in studies) {
    file <- file.path(results_folder, sprintf("study_%s.rds", study))
    problems <- c(problems, check_results_file(file, study, flatcor_commit, smoke))
  }
  problems <- c(problems, sprintf("figure script %s stopped", failed_figures))

  if (length(problems) == 0) {
    message("Check of the outputs: every results file and figure is there, ",
            "every number is a number.")
    return(invisible(TRUE))
  }
  report <- paste0("Check of the outputs found ", length(problems), " problem(s):\n",
                   paste0("  - ", problems, collapse = "\n"))
  if (smoke) stop(report, call. = FALSE)
  message(report)
  invisible(FALSE)
}
