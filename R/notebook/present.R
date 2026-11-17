## =====================================================================
##  Presenting the studies: the helpers of validation.ipynb
##
##  In one sentence: run the studies (optionally), read each study's
##  results/studies/study_NN.rds, and show its checks as a PASS/FAIL
##  table, its main tables formatted with gt (every estimate with its Monte
##  Carlo SE in brackets), and its figures.
##
##  Nothing here computes a result: the study scripts do that.
## =====================================================================
library(gt)

## Jupyter may start R in an ASCII locale; the checks contain names such as
## García-Donato, so use UTF-8 here and in the study processes started from here
if (!isTRUE(l10n_info()$`UTF-8`))
  for (locale in c("C.UTF-8", "en_US.UTF-8", "en_GB.UTF-8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break

## The four studies, in the order of docs/appendix/C_study_design.md
study_files <- c("01" = "01_prior_sampling", "02" = "02_bayes_factors", "03" = "03_calibration",
                 "04" = "04_default_prior_criteria", "2b" = "2b_estimation")
study_questions <- c("01" = "Does the code draw the prior it claims?",
                     "02" = "Are the Bayes factors computed correctly?",
                     "03" = "Are the posterior and the Bayes factor calibrated?",
                     "04" = "Does each prior meet the default-prior criteria?",
                     "2b" = "How close are the posterior estimates to a known truth?")

## Where the results are: the tracked folder, or cache/ for the smoke test
## (the same rule as R/shared/settings.R)
results_folder_of <- function(smoke) if (smoke) file.path("cache", "smoke_results") else "results"
cache_folder_of   <- function(study, smoke) file.path("cache", sprintf("study_%s%s", study, if (smoke) "_smoke" else ""))

## ---- Running ------------------------------------------------------------------------
## Runs each study in its own R process, as run_all.R does; with
## fresh = TRUE its cache folder is deleted first. Returns one row per study.
run_studies <- function(studies, fresh = FALSE, smoke = FALSE, cores = NA) {
  dir.create(file.path("cache", "logs"), recursive = TRUE, showWarnings = FALSE)
  environment <- c(paste0("VALIDATION_SMOKE=", if (smoke) "1" else "0"),
                   paste0("LC_CTYPE=", Sys.getlocale("LC_CTYPE")),
                   if (!is.na(cores)) paste0("VALIDATION_CORES=", cores))
  do.call(rbind, lapply(studies, function(study) {
    if (fresh) unlink(cache_folder_of(study, smoke), recursive = TRUE)
    log_file <- file.path("cache", "logs", paste0(study_files[[study]], ".log"))
    message(sprintf("[%s] study %s ...", format(Sys.time(), "%H:%M"), study))
    started <- Sys.time()
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      file.path("R", "studies", paste0(study_files[[study]], ".R")),
                      stdout = log_file, stderr = log_file, env = environment)
    data.frame(study = study, exit_status = status,
               minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
               log = log_file)
  }))
}

## ---- Reading ------------------------------------------------------------------------
read_study <- function(study, smoke = FALSE) {
  file <- file.path(results_folder_of(smoke), "studies", sprintf("study_%s.rds", study))
  if (file.exists(file)) readRDS(file) else NULL
}

## One row per study: its question, its checks and the run it comes from
overview_of <- function(studies, smoke = FALSE) {
  do.call(rbind, lapply(studies, function(study) {
    result <- read_study(study, smoke)
    if (is.null(result))
      return(data.frame(study = study, question = study_questions[[study]], passed = NA_integer_,
                        failed = NA_integer_, finished = "not run", commit = "", machine = "",
                        cores = NA_integer_, minutes = NA_real_))
    data.frame(study = study, question = study_questions[[study]],
               passed = sum(result$checks$passed), failed = sum(!result$checks$passed),
               finished = result$setup$finished, commit = result$setup$commit,
               machine = result$setup$processor %||% result$setup$machine %||% "",
               cores = as.integer(result$setup$cores),
               minutes = as.numeric(result$setup$minutes %||% NA))
  }))
}

## ---- Showing ------------------------------------------------------------------------
## gt tables are shown as HTML in Jupyter, and printed as text elsewhere
display_gt <- function(table) {
  if (isTRUE(getOption("jupyter.in_kernel")))
    IRdisplay::display_html(gt::as_raw_html(table, inline_css = TRUE))
  else print(table)
  invisible(table)
}

## Names in the data as they are written in the text
prior_labels <- c(flat_partial_r = "flat partial r", flat_r = "flat R", flat_r_squared = "flat R-squared")
column_labels <- c(k = "*k*", n = "*n*", q = "*q*", p_value = "*p*", mcse = "MCSE",
                   n_sim = "*n*<sub>sim</sub>", chi_square = "*χ*²", z = "*z*",
                   largest_abs_z = "Largest |*z*|", est = "Estimate")
abbreviations <- c(se = "SE", mcse = "MCSE", sd = "SD", mc = "MC", ks = "KS", ess = "ESS",
                   bf = "BF", jzs = "JZS")
label_of <- function(name) {
  if (name %in% names(column_labels)) return(column_labels[[name]])
  words <- strsplit(name, "_", fixed = TRUE)[[1]]
  known <- words %in% names(abbreviations)
  words[known] <- abbreviations[words[known]]
  label <- paste(words, collapse = " ")
  paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
}

## The Monte Carlo SE column that belongs to an estimate column, by the
## naming of the study scripts: "mcse" after its estimate; "mcse_<x>" with
## "<x>" or "..._<x>"; "<x>_mcse" with "<x>", "<x>_mean" or "<x>_...".
mcse_pairs <- function(columns) {
  pairs <- list()
  for (position in seq_along(columns)) {
    name <- columns[position]
    estimate <- if (name == "mcse") columns[position - 1] else
      if (startsWith(name, "mcse_")) {
        x <- sub("^mcse_", "", name)
        c(columns[columns == x], columns[endsWith(columns, paste0("_", x)) & columns != name])[1]
      } else if (endsWith(name, "_mcse")) {
        x <- sub("_mcse$", "", name)
        c(columns[columns %in% c(x, paste0(x, "_mean"))],
          columns[startsWith(columns, paste0(x, "_")) & columns != name])[1]
      } else NA
    if (length(estimate) == 1 && !is.na(estimate) && !startsWith(estimate, "mcse"))
      pairs[[estimate]] <- name
  }
  pairs
}

## Decimals for a column of numbers, from their typical size
decimals_for <- function(values) {
  typical <- stats::median(abs(values[is.finite(values) & values != 0]))
  if (!is.finite(typical)) 2 else if (typical >= 100) 1 else if (typical >= 1) 2 else 3
}

## A data frame as a gt table: priors and columns named as in the text, whole
## numbers as integers, every estimate with its Monte Carlo SE in brackets,
## and a scrolling box for long tables
show_table <- function(data, title = NULL, note = NULL) {
  data <- as.data.frame(data)
  for (column in names(data)) {
    if (is.character(data[[column]]) || is.factor(data[[column]])) {
      values <- as.character(data[[column]])
      known <- values %in% names(prior_labels)
      values[known] <- prior_labels[values[known]]
      data[[column]] <- values
    }
    if (is.list(data[[column]])) data[[column]] <- NULL
  }
  pairs <- mcse_pairs(names(data))
  table <- gt(data)
  for (column in names(data)) {
    values <- data[[column]]
    if (!is.numeric(values) || column %in% unlist(pairs)) next
    finite <- values[is.finite(values)]
    if (length(finite) > 0 && all(finite == round(finite)) && !column %in% names(pairs)) {
      table <- fmt_integer(table, columns = all_of(column))
      next
    }
    decimals <- decimals_for(values)
    table <- fmt_number(table, columns = all_of(c(column, pairs[[column]])), decimals = decimals)
    if (!is.null(pairs[[column]]))
      table <- cols_merge(table, columns = all_of(c(column, pairs[[column]])), pattern = "{1} ({2})")
  }
  table <- cols_label(table, .list = setNames(lapply(names(data), function(n) md(label_of(n))), names(data)))
  table <- tab_options(table, table.font.size = px(13), data_row.padding = px(3),
                       heading.align = "left", table.align = "left")
  table <- sub_missing(table, missing_text = "–")
  if (!is.null(title)) table <- tab_header(table, title = title)
  notes <- c(if (length(pairs) > 0) "Monte Carlo standard errors in brackets.", note)
  if (length(notes) > 0) table <- tab_source_note(table, paste(notes, collapse = " "))
  if (nrow(data) > 20) table <- tab_options(table, container.height = px(420), container.overflow.y = TRUE)
  display_gt(table)
}

## The checks of a study, PASS in green and FAIL in red
show_checks <- function(result) {
  checks <- result$checks
  if (nrow(checks) == 0)
    return(show_text("No pass/fail checks: each estimate is shown with its reference and its Monte Carlo SE."))
  data <- data.frame(result = ifelse(checks$passed, "PASS", "FAIL"), check = checks$check,
                     details = checks$details)
  table <- gt(data) |>
    cols_label(result = "", check = "Check", details = "Details") |>
    tab_style(cell_text(color = "#1a7f37", weight = "bold"),
              cells_body(columns = result, rows = result == "PASS")) |>
    tab_style(cell_text(color = "#cf222e", weight = "bold"),
              cells_body(columns = result, rows = result == "FAIL")) |>
    tab_header(title = sprintf("Checks: %d of %d pass", sum(checks$passed), nrow(checks))) |>
    tab_options(table.font.size = px(13), data_row.padding = px(3), heading.align = "left",
                table.align = "left")
  if (isTRUE(result$setup$smoke))
    table <- tab_source_note(table, "Smoke test: with 10 datasets per cell the checks mean nothing.")
  display_gt(table)
}

## A line of text in its place among the tables (cat() output can arrive after them in Jupyter)
show_text <- function(text)
  if (isTRUE(getOption("jupyter.in_kernel"))) IRdisplay::display_markdown(text) else cat(text, "\n")

show_figures <- function(result) {
  for (file in result$figures) {
    if (!file.exists(file)) { message("Missing figure: ", file); next }
    if (isTRUE(getOption("jupyter.in_kernel"))) IRdisplay::display_png(file = file)
    else message("Figure: ", file)
  }
}

## Everything of one study: the run it comes from, its checks, the chosen
## tables and its figures. `tables` maps a title to the name of a table in
## the results file, or to a function of all its tables that returns one.
show_study <- function(study, tables = list(), smoke = FALSE) {
  result <- read_study(study, smoke)
  if (is.null(result)) {
    show_text(sprintf("Study %s has no results yet (`%s/studies/study_%s.rds`). Run it, or set `mode` to `\"cached\"`.",
                      study, results_folder_of(smoke), study))
    return(invisible(NULL))
  }
  show_text(sprintf("Run of %s, commit %s, flatcor %s, %s cores%s%s%s.", result$setup$finished,
                    result$setup$commit, result$setup$flatcor, result$setup$cores,
                    if (is.na(result$setup$processor %||% NA)) "" else paste0(" of ", result$setup$processor),
                    if (is.null(result$setup$minutes)) "" else sprintf(", %.0f minutes", result$setup$minutes),
                    if (isTRUE(result$setup$smoke)) ", **smoke test**" else ""))
  show_checks(result)
  for (title in names(tables)) {
    data <- if (is.function(tables[[title]])) tables[[title]](result$tables) else result$tables[[tables[[title]]]]
    if (is.null(data)) { message("No table '", tables[[title]], "' in study ", study); next }
    show_table(data, title)
  }
  show_figures(result)
  invisible(result)
}
