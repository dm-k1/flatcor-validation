## =====================================================================
##  The static check: seconds, nothing sampled, nothing run but the
##  definitions.
##
##    Rscript R/shared/static_check.R
##
##    1. every R file parses
##    2. the pipeline is valid (targets::tar_validate()) and its steps can be
##       listed (targets::tar_manifest(), which runs _targets.R and with it
##       every design.R)
##    3. every function the pipeline's code calls exists: in its own files,
##       in a package it loads, or as a package::function that package has
##    4. the code keeps to the declared way (section 4 below): every
##       expensive call inside a seam the dry run replays, every file a step needs or makes declared
##       as a step, every results file read made by a step
##
##  It stops with an error listing the problems. The dry run (R/dry_run/replay.R)
##  tests what this cannot: the tables handed from step to step.
## =====================================================================

found <- new.env()
found$problems <- character(0)
note <- function(...) found$problems <- unique(c(found$problems, sprintf(...)))

## ---- 1. Every R file parses --------------------------------------------------------
r_files <- c(
  "_targets.R", "run_all.R",
  list.files("R", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
)
for (file in r_files) {
  parsed <- tryCatch(parse(file), error = function(e) e)
  if (inherits(parsed, "error")) note("%s does not parse: %s", file, conditionMessage(parsed))
}
cat(sprintf("1. %d R files parsed\n", length(r_files)))

## ---- 2. The pipeline --------------------------------------------------------------
manifest <- tryCatch(
  {
    targets::tar_validate()
    targets::tar_manifest(fields = c("name", "command"))
  },
  error = function(e) {
    note("the pipeline is not valid: %s", conditionMessage(e))
    NULL
  }
)
if (!is.null(manifest)) cat(sprintf("2. the pipeline is valid: %d steps\n", nrow(manifest)))

## ---- 3. Every function called exists ------------------------------------------------
## The pipeline's code as the workers have it: the files _targets.R sources,
## with the packages it loads
studies <- basename(list.dirs(file.path("R", "2_studies"), recursive = FALSE))
pipeline_files <- c(
  file.path("R", "shared", c(
    "random_streams.R", "pipeline.R", "expectations.R", "results.R", "per_dataset.R",
    "comparison_methods.R", "rounds.R", "check_results.R", "workers.R", "queue.R"
  )),
  list.files(file.path("R", "1_data"), pattern = "[.]R$", full.names = TRUE),
  file.path("R", "2_studies", studies, "design.R"),
  file.path("R", "2_studies", studies, "code.R")
)
# packages the pipeline loads (_targets.R), and those only the figures load,
# whose helpers live in comparison_methods.R
for (package in c("targets", "flatcor", "dplyr", "ggplot2", "grid")) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}
code <- new.env()
for (file in pipeline_files) sys.source(file, envir = code)

## Call heads in an expression: the names called (f in f(...)), and pkg::f
## pairs, leaving out names that are defined in it (arguments and local
## assignments) so local helpers count as defined
calls_in <- function(expression) {
  called <- character(0)
  namespaced <- list()
  defined <- character(0)
  walk <- function(e) {
    if (is.function(e)) {
      defined <<- c(defined, names(formals(e)))
      walk(body(e))
      return()
    }
    if (!is.call(e)) return()
    head <- e[[1]]
    if (is.symbol(head)) {
      name <- as.character(head)
      if (name %in% c("<-", "=", "<<-", "assign") && length(e) >= 2) {
        target <- e[[2]]
        if (is.symbol(target)) defined <<- c(defined, as.character(target))
        if (is.character(target)) defined <<- c(defined, target)
      }
      if (name == "function") defined <<- c(defined, names(e[[2]]))
      if (name == "for") defined <<- c(defined, as.character(e[[2]]))
      # inside a formula, names are columns, not calls
      if (name == "~") return()
      called <<- c(called, name)
    } else if (is.call(head) && deparse(head[[1]]) %in% c("::", ":::")) {
      namespaced[[length(namespaced) + 1]] <<- c(as.character(head[[2]]), as.character(head[[3]]))
    } else {
      walk(head)
    }
    for (argument in as.list(e)[-1]) if (!missing(argument)) walk(argument)
  }
  walk(expression)
  list(called = setdiff(unique(called), defined), namespaced = unique(namespaced))
}

## Functions that work on columns inside dplyr verbs, called by name there
column_helpers <- c("n", "across", "everything", "where", "starts_with", "ends_with")

check_calls <- function(expression, where) {
  calls <- calls_in(expression)
  for (name in setdiff(calls$called, column_helpers)) {
    if (!exists(name, envir = code, mode = "function")) {
      note("%s calls %s(), which is not defined", where, name)
    }
  }
  for (pair in calls$namespaced) {
    if (!requireNamespace(pair[1], quietly = TRUE)) {
      note("%s calls %s::%s, but %s is not installed", where, pair[1], pair[2], pair[1])
    } else if (!exists(pair[2], envir = asNamespace(pair[1]))) {
      note("%s calls %s::%s, which %s does not have", where, pair[1], pair[2], pair[1])
    }
  }
}

functions <- Filter(function(name) is.function(code[[name]]), ls(code))
for (name in functions) check_calls(code[[name]], sprintf("%s()", name))
if (!is.null(manifest)) {
  for (i in seq_len(nrow(manifest))) {
    check_calls(str2lang(manifest$command[[i]]), sprintf("step %s", manifest$name[[i]]))
  }
}
cat(sprintf("3. calls checked in %d functions and the pipeline's steps\n", length(functions)))

## ---- 4. The declared way -------------------------------------------------------------
## The pipeline's code keeps to three rules, so that targets knows what every
## step needs and makes, and the dry run can stand in for every expensive
## call:
##   a. every expensive call (Stan, bridge sampling, SBC) is inside one of
##      the seams the dry run replays (replayed_functions in
##      R/dry_run/replay.R), and every seam listed there exists
##   b. files are read and written only where it is declared: a study's results
##      file (save_study_results_file(), a "file" step), the plans of the
##      rounds and the pieces the queue runs (both in the store), the
##      recordings, the compile lock and the final check. Any other file a
##      step needs or makes is a step of its own (tar_target(..., format = "file"))
##   c. every results file a figure script or run_all.R reads is made by a
##      step: results_<study>, format = "file"; the figure scripts read them
##      only through read_results()

expensive_calls <- c(
  "bayes_factor_models", "fit_regression", "compute_SBC", "bridge_sampler",
  "cmdstan_model", "sample" # sample: model$sample(), a CmdStan fit
)
file_calls <- c(
  "readRDS", "saveRDS", "load", "save", "read.csv", "write.csv", "readLines", "writeLines",
  "read_csv", "write_csv", "fread", "fwrite", "file.exists", "file.rename", "file.copy",
  "unlink", "dir.create", "list.files", "ggsave", "png", "pdf"
)
allowed_file_use <- c(
  "save_study_results_file", "read_plan", "write_plan",
  "run_piece", "read_piece", "run_rounds", "as_read_back", "code_files",
  "compile_flatcor_models", "check_results_file", "check_results",
  "write_per_dataset", "read_per_dataset"
)

## The seams, as the dry run's file lists them (sourced on its own, it
## replaces nothing)
overlay <- new.env()
sys.source(file.path("R", "dry_run", "replay.R"), envir = overlay)
seams <- overlay$replayed_functions
for (skipped in overlay$skipped_in_dry_run) {
  if (!exists(skipped, envir = code, mode = "function")) {
    note("the dry run skips %s(), which no longer exists (R/dry_run/replay.R)", skipped)
  }
}
# flatcor's fit is replaced by name: called as fit_regression(), never as
# flatcor::fit_regression(), which would pass the replacement by
flatcor_seams <- intersect(seams, getNamespaceExports("flatcor"))
for (seam in setdiff(seams, flatcor_seams)) {
  if (!exists(seam, envir = code, mode = "function")) {
    note("the dry run replays %s(), which no longer exists (R/dry_run/replay.R)", seam)
  }
}

## The calls in an expression as (name, how) pairs: name is f for f(...),
## pkg::f(...) and x$f(...); how is "plain", "::" or "$"
calls_with_context <- function(expression) {
  found <- list()
  walk <- function(e) {
    if (is.function(e)) return(walk(body(e)))
    if (!is.call(e)) return()
    head <- e[[1]]
    if (is.symbol(head)) {
      found[[length(found) + 1]] <<- c(name = as.character(head), how = "plain")
    } else if (is.call(head) && deparse(head[[1]]) %in% c("::", ":::", "$")) {
      how <- if (deparse(head[[1]]) == "$") "$" else "::"
      found[[length(found) + 1]] <<- c(name = as.character(head[[3]]), how = how)
    } else {
      walk(head)
    }
    for (argument in as.list(e)[-1]) if (!missing(argument)) walk(argument)
  }
  walk(expression)
  found
}

for (name in functions) {
  calls <- calls_with_context(code[[name]])
  for (call in calls) {
    called <- call[["name"]]
    how <- call[["how"]]
    # the base sample() also exists: only x$sample() is Stan's
    stan_sample <- called == "sample" && how == "$"
    replaced_by_name <- called %in% flatcor_seams && how == "plain"
    is_expensive <- called %in% setdiff(expensive_calls, "sample") || stan_sample
    if (is_expensive && !(name %in% seams) && !replaced_by_name) {
      note(
        paste(
          "%s() calls %s() outside the seams the dry run replays: move it into one,",
          "or add %s() to replayed_functions with its replacement (R/dry_run/replay.R)"
        ),
        name, called, name
      )
    }
    if (called %in% file_calls && !(name %in% allowed_file_use)) {
      note(
        "%s() uses %s(): a file a step needs or makes is a step of its own (format = \"file\")",
        name, called
      )
    }
  }
}

if (!is.null(manifest)) {
  formats <- targets::tar_manifest(fields = c("name", "format"))
  is_results_file_step <- function(study) {
    step <- formats[formats$name == paste0("results_", study), ]
    nrow(step) == 1 && identical(step$format[[1]], "file")
  }
  figure_scripts <- list.files(file.path("R", "3_figures"), pattern = "[.]R$", full.names = TRUE)
  for (script in figure_scripts) {
    for (call in calls_with_context(as.call(c(as.name("{"), as.list(parse(script)))))) {
      if (call[["name"]] %in% c("readRDS", "load", "read.csv", "fread")) {
        note("%s reads a file with %s(): read results with read_results()", script, call[["name"]])
      }
    }
    studies_read <- regmatches(
      readLines(script),
      regexpr("read_results[(]\"[^\"]+\"[)]", readLines(script))
    )
    for (study in sub("read_results[(]\"([^\"]+)\"[)]", "\\1", studies_read)) {
      if (!is_results_file_step(study)) {
        note(
          "%s reads study %s, but no step results_%s (format = \"file\") makes it",
          script, study, study
        )
      }
    }
  }
  # run_all.R reads, at the end, the results file of every study that runs in rounds
  run_all_lines <- grep("study = \"[0-9][a-z]\"", readLines("run_all.R"), value = TRUE)
  run_all_studies <- unique(sub(".*study = \"([0-9][a-z])\".*", "\\1", run_all_lines))
  for (study in run_all_studies) {
    if (!is_results_file_step(study)) {
      note(
        "run_all.R reads study %s, but no step results_%s (format = \"file\") makes it",
        study, study
      )
    }
  }
}
cat("4. the declared way: expensive calls in seams, files declared, results made by steps\n")

## ---- The verdict ----------------------------------------------------------------------
if (length(found$problems) > 0) {
  stop(
    sprintf("The static check found %d problem(s):\n", length(found$problems)),
    paste0("  - ", found$problems, collapse = "\n"),
    call. = FALSE
  )
}
cat("The static check found no problems.\n")
