## =====================================================================
##  The dry run: the expensive calls recorded once, then replayed
##
##  The dry run (VALIDATION_DRY=1 Rscript run_all.R) tests the plumbing of
##  the whole validation in a minute or two: every step, every hand-over of
##  a table from one step to the next, the rounds, every results file, every
##  figure and the check at the end. It runs the smoke test's design, in one
##  R process, with the calls to Stan and bridge sampling replaced by what
##  those same calls returned in a real smoke test.
##
##    VALIDATION_SMOKE=1 VALIDATION_RECORD=1 Rscript run_all.R
##        a smoke test that also records what each expensive call returned
##        (cache/dry_recordings/); once, and again when flatcor or what an
##        expensive call returns changes
##    VALIDATION_DRY=1 Rscript run_all.R
##        the dry run, from those recordings
##
##  This file changes no study's code. It is sourced after the studies'
##  files, and only in a recording run or the dry run (_targets.R, run_all.R);
##  a normal run never sees it. It replaces a few functions by versions that
##  record or replay the original: the seams, listed in replayed_functions.
##  Every expensive call in the studies' code must be inside one of them
##  (the static check enforces it: R/shared/static_check.R).
##
##  What the dry run tests: everything around the seams. What it does not:
##  the code inside them (Stan, bridge sampling, sample_setting_1a()'s
##  lengthening of the chains, ...), and whether the numbers are right; the
##  smoke test does that.
##
##  A recording is found by its call's exact inputs. A call with inputs not
##  seen in the recording run (a changed design or data step) gets a
##  recording of the same shape (same prior, numbers of predictors, ...), so
##  the dry run still runs, and says so once. With no recording of that
##  shape it stops and says what to do.
## =====================================================================

dry <- identical(Sys.getenv("VALIDATION_DRY"), "1")
recording <- identical(Sys.getenv("VALIDATION_RECORD"), "1")
recordings_folder <- file.path("cache", "dry_recordings")

## The seams: the functions this file replaces. fit_regression is flatcor's,
## called by the studies without flatcor::, so a function of that name here
## comes first.
replayed_functions <- c(
  "log_bayes_factor_between", # comparison_methods.R: flatcor's Bayes factors (1b, 3a, 3b)
  "fit_regression",           # flatcor's fit (2b)
  "sample_setting_1a",        # 1a: the prior draws of one setting
  "sbc_batch_2a"              # 2a: the SBC fits of one batch
)

## Skipped in the dry run, which samples nothing: compiling flatcor's models
## (R/shared/pipeline.R). The static check fails when one of these no longer
## exists, so a rename never silently brings compiling back into the dry run.
skipped_in_dry_run <- c("compile_flatcor_models", "flatcor_stan_model")


## ---- Record or replay one call --------------------------------------------------------

## The value of `expression`, run and saved (recording) or read from a
## recording (dry run).
##   name   the kind of call, a folder in cache/dry_recordings/
##   shape  a named list of what fixes the shape of the value (prior, k, ...)
##   inputs everything the value depends on (data, seed, settings), hashed
##          to find the recording of this exact call; no model objects
replayable <- function(name, shape, inputs, expression) {
  if (!dry && !recording) {
    return(expression)
  }
  shape_key <- paste(
    c(name, sprintf("%s-%s", names(shape), vapply(shape, paste, "", collapse = "x"))),
    collapse = "_"
  )
  shape_key <- gsub("[^A-Za-z0-9_.-]", "-", shape_key)
  exact <- substr(digest::digest(inputs), 1, 16)
  folder <- file.path(recordings_folder, name)
  file <- file.path(folder, sprintf("%s__%s.rds", shape_key, exact))

  if (recording) {
    value <- expression
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    # written under another name, then renamed: a worker that stops halfway
    # leaves no broken recording
    partial <- tempfile(tmpdir = folder, fileext = ".partial")
    saveRDS(value, partial)
    file.rename(partial, file)
    return(value)
  }

  # the dry run
  if (file.exists(file)) {
    return(readRDS(file))
  }
  recorded <- list.files(folder, pattern = "[.]rds$")
  same_shape <- recorded[startsWith(recorded, paste0(shape_key, "__"))]
  if (length(same_shape) == 0) {
    stop(
      "Dry run: no recording of ", shape_key, " (", folder, "). Record them with a ",
      "smoke test: VALIDATION_SMOKE=1 VALIDATION_RECORD=1 Rscript run_all.R",
      call. = FALSE
    )
  }
  if (!isTRUE(getOption(paste0("validation_replay_warned_", name)))) {
    message(
      "Dry run: ", name, " called with inputs the recording run did not have; ",
      "using recordings of the same shape (the design or a data step changed ",
      "since the recordings were made)"
    )
    options(stats::setNames(list(TRUE), paste0("validation_replay_warned_", name)))
  }
  # the same recording for the same inputs, different ones for different inputs
  pick <- strtoi(substr(exact, 1, 7), 16L) %% length(same_shape) + 1
  readRDS(file.path(folder, same_shape[[pick]]))
}


## ---- The seams --------------------------------------------------------------------------

## Where this file is sourced: the pipeline's functions are there
overlay_env <- environment()

## Replaces function `name` (if it is defined here, and not replaced already)
## by make(original)
replace_function <- function(name, make) {
  if (!exists(name, envir = overlay_env, inherits = FALSE)) return(invisible())
  original <- get(name, envir = overlay_env)
  if (isTRUE(attr(original, "replayed"))) return(invisible())
  replacement <- make(original)
  attr(replacement, "replayed") <- TRUE
  assign(name, replacement, envir = overlay_env)
}

## flatcor's Bayes factor between two models (the three priors; the
## references are cheap and run as they are)
replace_function("log_bayes_factor_between", function(original) {
  function(method, predictors_a, predictors_b, outcome, seed) {
    if (!(method %in% names(contenders))) {
      return(original(method, predictors_a, predictors_b, outcome, seed))
    }
    replayable(
      "bayes_factor",
      shape = list(
        prior = contenders[[method]], k_a = ncol(predictors_a), k_b = ncol(predictors_b)
      ),
      inputs = list(method, predictors_a, predictors_b, outcome, seed, bayes_factor_draws),
      original(method, predictors_a, predictors_b, outcome, seed)
    )
  }
})

## 1a: the prior draws of one setting (`model` is not part of the inputs:
## in the dry run there is none)
replace_function("sample_setting_1a", function(original) {
  function(setting, predictors, model) {
    replayable(
      "sample_setting_1a",
      shape = list(prior = setting$prior, design = setting$design),
      inputs = list(setting, predictors, start_draws_1a),
      original(setting, predictors, model)
    )
  }
})

## 2a: the SBC fits of one batch
replace_function("sbc_batch_2a", function(original) {
  function(setting, replications, model) {
    replayable(
      "sbc_batch_2a",
      shape = list(prior = setting$prior, k = setting$k, fits = length(replications)),
      inputs = list(setting, replications, iterations_2a),
      original(setting, replications, model)
    )
  }
})

## flatcor's fit (2b). What is kept of the fit: the parts the studies use,
## its draws and its divergent transitions.
if (exists("estimates_2b", envir = overlay_env, inherits = FALSE)) {
  fit_regression <- function(formula, data, prior, ..., seed) {
    replayable(
      "fit_regression",
      shape = list(prior = prior, k = length(all.vars(formula)) - 1),
      inputs = list(deparse(formula), data, prior, list(...), seed),
      {
        fit <- flatcor::fit_regression(formula, data, prior = prior, ..., seed = seed)
        list(draws = fit$draws, divergent = fit$divergent)
      }
    )
  }
  attr(fit_regression, "replayed") <- TRUE
}

## The dry run samples nothing, so compiles no model
if (dry) {
  replace_function("compile_flatcor_models", function(original) {
    function() invisible(TRUE)
  })
  replace_function("flatcor_stan_model", function(original) {
    function() "dry run: no model"
  })
}
