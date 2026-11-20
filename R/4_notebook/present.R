## =====================================================================
##  Presenting the studies: the helpers of validation.ipynb
##
##  In one sentence: read each study's results/studies/study_<cell>.rds and
##  show its main tables formatted with gt (every estimate against what is
##  expected, with its Monte Carlo SE in brackets), any tripwire that failed,
##  and its figures (drawn with make_figures.R when they are not there yet).
##
##  Nothing here runs a study or computes a result: run_all.R runs the
##  studies (on a big machine), the study scripts compute.
## =====================================================================
library(gt)

## Jupyter may start R in an ASCII locale; the results contain names such as
## García-Donato, so use UTF-8 here and in the figure scripts started from here
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (locale in c("C.UTF-8", "en_US.UTF-8", "en_GB.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break
  }
}

## The six studies, one per cell of the table in docs/working/REDESIGN.md:
## heading (1 prior sampling, 2 posterior estimation, 3 Bayes factors) and row
## (a: computed correctly, b: behaves like a default)
study_files <- c(
  "1a" = "1a_prior_sampling", "1b" = "1b_invariance", "2a" = "2a_posterior",
  "2b" = "2b_estimation", "3a" = "3a_bayes_factors", "3b" = "3b_default_criteria"
)
study_questions <- c(
  "1a" = "Does the software draw the prior it claims?",
  "1b" = "Is the Bayes factor invariant to the units and coding of the predictors?",
  "2a" = "Is the posterior computed correctly (simulation-based calibration)?",
  "2b" = "How close are the posterior estimates to a known truth?",
  "3a" = "Are the Bayes factors right (calibration; bridge sampling's error)?",
  "3b" = "Do the Bayes factors meet the default-prior criteria?"
)

## Where the results are: the tracked folder, or cache/ for the smoke test
## (the same rule as R/shared/results.R)
results_folder_of <- function(smoke) if (smoke) file.path("cache", "smoke_results") else "results"

## ---- Reading ------------------------------------------------------------------------
read_study <- function(study, smoke = FALSE) {
  file <- file.path(results_folder_of(smoke), "studies", sprintf("study_%s.rds", study))
  if (file.exists(file)) readRDS(file) else NULL
}

## One row per study: its question and the run it comes from
overview_of <- function(studies, smoke = FALSE) {
  do.call(rbind, lapply(studies, function(study) {
    result <- read_study(study, smoke)
    if (is.null(result)) {
      return(data.frame(
        study = study, question = study_questions[[study]], finished = "not run",
        commit = "", flatcor = "", machine = ""
      ))
    }
    data.frame(
      study = study, question = study_questions[[study]],
      finished = result$setup$finished, commit = result$setup$commit,
      flatcor = result$setup$flatcor, machine = result$setup$machine %||% ""
    )
  }))
}

## ---- Showing ------------------------------------------------------------------------
## gt tables are shown as HTML in Jupyter, and printed as text elsewhere
display_gt <- function(table) {
  if (isTRUE(getOption("jupyter.in_kernel"))) {
    IRdisplay::display_html(gt::as_raw_html(table, inline_css = TRUE))
  } else {
    print(table)
  }
  invisible(table)
}

## Names in the data as they are written in the text
prior_labels <- c(
  flat_partial_r = "flat partial r",
  flat_r = "flat R",
  flat_r_squared = "flat R-squared"
)
column_labels <- c(
  k = "*k*", n = "*n*", q = "*q*", p_value = "*p*", mcse = "MCSE",
  n_sim = "*n*<sub>sim</sub>", chi_square = "*χ*²", z = "*z*",
  largest_abs_z = "Largest |*z*|", est = "Estimate"
)
abbreviations <- c(
  se = "SE", mcse = "MCSE", sd = "SD", mc = "MC", ks = "KS", ess = "ESS",
  bf = "BF", jzs = "JZS"
)
label_of <- function(name) {
  if (name %in% names(column_labels)) {
    return(column_labels[[name]])
  }
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
    estimate <- if (name == "mcse") {
      columns[position - 1]
    } else if (startsWith(name, "mcse_")) {
      x <- sub("^mcse_", "", name)
      c(columns[columns == x], columns[endsWith(columns, paste0("_", x)) & columns != name])[1]
    } else if (endsWith(name, "_mcse")) {
      x <- sub("_mcse$", "", name)
      c(
        columns[columns %in% c(x, paste0(x, "_mean"))],
        columns[startsWith(columns, paste0(x, "_")) & columns != name]
      )[1]
    } else {
      NA
    }
    if (length(estimate) == 1 && !is.na(estimate) && !startsWith(estimate, "mcse")) {
      pairs[[estimate]] <- name
    }
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
    if (!is.null(pairs[[column]])) {
      table <- cols_merge(
        table,
        columns = all_of(c(column, pairs[[column]])),
        pattern = "{1} ({2})"
      )
    }
  }
  labels <- setNames(lapply(names(data), function(n) md(label_of(n))), names(data))
  table <- cols_label(table, .list = labels)
  table <- tab_options(table,
    table.font.size = px(13), data_row.padding = px(3),
    heading.align = "left", table.align = "left"
  )
  table <- sub_missing(table, missing_text = "–")
  if (!is.null(title)) table <- tab_header(table, title = title)
  notes <- c(if (length(pairs) > 0) "Monte Carlo standard errors in brackets.", note)
  if (length(notes) > 0) table <- tab_source_note(table, paste(notes, collapse = " "))
  if (nrow(data) > 20) {
    table <- tab_options(table, container.height = px(420), container.overflow.y = TRUE)
  }
  display_gt(table)
}

## A line of text in its place among the tables (cat() output can arrive after them in Jupyter)
show_text <- function(text) {
  if (isTRUE(getOption("jupyter.in_kernel"))) IRdisplay::display_markdown(text) else cat(text, "\n")
}

## The figures of a study: the PNGs its script in R/3_figures/ wrote (named
## after the study, as results/figures/3a_*.png). When there are none yet
## (results downloaded from a run, say), make_figures.R draws them first.
show_figures <- function(result) {
  smoke <- isTRUE(result$setup$smoke)
  folder <- if (smoke) file.path("cache", "smoke_figures") else file.path("results", "figures")
  pattern <- paste0("^", result$study, "_.*[.]png$")
  has_script <- any(startsWith(list.files(file.path("R", "3_figures")), paste0(result$study, "_")))
  if (has_script && length(list.files(folder, pattern = pattern)) == 0) {
    status <- system2(file.path(R.home("bin"), "Rscript"), c("make_figures.R", result$study),
      env = c(
        paste0("VALIDATION_SMOKE=", if (smoke) "1" else "0"),
        paste0("LC_CTYPE=", Sys.getlocale("LC_CTYPE"))
      )
    )
    if (status != 0) message("make_figures.R ", result$study, " stopped (see above)")
  }
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  for (file in files) {
    if (isTRUE(getOption("jupyter.in_kernel"))) {
      IRdisplay::display_png(file = file)
    } else {
      message("Figure: ", file)
    }
  }
}

## Divergent transitions in one sentence, from a study's divergences table:
## `total` and `with` name its columns of fits (or Bayes factors) in all and
## of those with divergent transitions; `what` says what they are
divergence_sentence <- function(table, total, with, what = "fits") {
  table <- as.data.frame(table)
  sprintf(
    "Divergent transitions: in %s of %s %s (%s transitions in all).",
    format(sum(table[[with]]), big.mark = ","), format(sum(table[[total]]), big.mark = ","),
    what, format(sum(table$divergent_transitions), big.mark = ",")
  )
}

## 2b in short: per method and quantity, the mean and the range of the
## estimates over the cells (the full table is under Details)
summary_2b <- function(tables) {
  as.data.frame(tables$estimates) |>
    dplyr::summarise(
      cells = dplyr::n(), mean = mean(estimate), lowest = min(estimate), highest = max(estimate),
      .by = c(quantity, method)
    ) |>
    dplyr::arrange(quantity, method)
}

## Everything of one study: the run it comes from, the chosen tables, a
## sentence or two and its figures. `tables` maps a title to the name of a
## table in the results file, or to a function of all its tables that
## returns one; `notes` are functions of all its tables that return a
## sentence (divergence_sentence(), say).
show_study <- function(study, tables = list(), notes = list(), smoke = FALSE) {
  result <- read_study(study, smoke)
  if (is.null(result)) {
    show_text(sprintf(
      paste(
        "Study %s has no results yet (`%s/studies/study_%s.rds`):",
        "run it with `Rscript run_all.R %s` (README, \"Running the validation\")."
      ),
      study, results_folder_of(smoke), study, study
    ))
    return(invisible(NULL))
  }
  show_text(sprintf(
    "Run of %s, commit %s, flatcor %s, %s%s. Each estimate is shown against what is expected, %s",
    result$setup$finished, result$setup$commit, result$setup$flatcor, result$setup$machine,
    if (isTRUE(result$setup$smoke)) ", **smoke test: its numbers mean nothing**" else "",
    "with its Monte Carlo SE (no pass/fail)."
  ))
  for (title in names(tables)) {
    data <- if (is.function(tables[[title]])) {
      tables[[title]](result$tables)
    } else {
      result$tables[[tables[[title]]]]
    }
    if (is.null(data)) {
      message("No table '", tables[[title]], "' in study ", study)
      next
    }
    show_table(data, title)
  }
  for (note in notes) show_text(note(result$tables))
  show_figures(result)
  invisible(result)
}

## The tables of a study that are not needed to see its result (under
## Details at the end of the notebook): tables as for show_study()
show_details <- function(study, tables, smoke = FALSE) {
  result <- read_study(study, smoke)
  if (is.null(result)) return(invisible(NULL))
  show_text(sprintf("**Study %s**", study))
  for (title in names(tables)) {
    data <- result$tables[[tables[[title]]]]
    if (!is.null(data)) show_table(data, title)
  }
  invisible(result)
}
