## =====================================================================
##  One build of flatcor's Bayes factor model, shared by every worker
##
##  flatcor builds its Stan model in two versions: for estimation, kept in
##  R's user cache folder and reused by every process, and with "model
##  methods" for the Bayes factors (bridge sampling evaluates the log density
##  from R). The second version flatcor builds again in every R process (about
##  45 seconds), as cmdstanr does: cmdstanr compiles the model methods with
##  Rcpp from a temporary file, so they are never reused. For daily use that
##  is right (once per session). In a simulation every worker pays it, every
##  run.
##
##  This file builds that version once, in cache/stan_shared_build/, and
##  gives every worker that build (a few seconds to load). It is the only
##  place that does so, and it reaches into the internals of two packages:
##    flatcor  its in-memory store of compiled models (compiled_models) and
##             the name it keeps the Bayes factor model under
##    cmdstanr how it builds the model methods (get_standalone_hpp,
##             get_cmdstan_flags, repair_path, model_methods.cpp) and where a
##             model keeps them (private field model_methods_env_)
##  renv pins both packages, so these do not change underneath a run. When a
##  check below finds one changed, the workers fall back to flatcor's own
##  build (slower, never wrong) and say so.
##
##  The explicit test (about two minutes the first time, a build per version):
##    Rscript R/shared/shared_stan_build.R
##  checks that every internal in shared_build_internals is still there
##  (naming similar names when one is gone: a rename?), then computes one Bayes factor
##  with flatcor's own build and one with the shared build, each in a fresh R
##  process, and compares them: they must be identical.
## =====================================================================

## Where the shared build is: one folder per flatcor, CmdStan and cmdstanr
## version (the model methods' C++ comes from cmdstanr), so an update builds
## anew instead of using an old build
shared_build_folder <- function() {
  file.path("cache", "stan_shared_build", paste(
    "flatcor", utils::packageVersion("flatcor"),
    "cmdstan", cmdstanr::cmdstan_version(),
    "cmdstanr", utils::packageVersion("cmdstanr"),
    sep = "_"
  ))
}

## The name flatcor keeps its Bayes factor model under (compile_stan_model())
shared_build_key <- "cmdstanr methods"


## ---- Is the shared build still possible? ----------------------------------------------

## Every internal the shared build relies on, declared once. The test walks
## this list, so a new reliance is one row here.
##   package  where it lives
##   name     its name there (a function, an object, a file, a field)
##   kind     function, environment, file (in the package's folders) or
##            private_field (of cmdstanr's CmdStanModel)
##   use      what the shared build does with it
shared_build_internals <- data.frame(
  package = c("flatcor", "flatcor", "cmdstanr", "cmdstanr", "cmdstanr", "cmdstanr", "cmdstanr"),
  name = c(
    "compiled_models", "compile_stan_model", "get_standalone_hpp", "get_cmdstan_flags",
    "repair_path", "include/model_methods.cpp", "model_methods_env_"
  ),
  kind = c(
    "environment", "function", "function", "function",
    "function", "file", "private_field"
  ),
  use = c(
    "where flatcor keeps its compiled models; the shared build is put there",
    "keeps the Bayes factor model under shared_build_key (checked in its code)",
    "the model's C++ code, for the model methods",
    "CmdStan's compiler flags, for the model methods",
    "the path to CmdStan's TBB library while compiling",
    "the C++ of the model methods",
    "where a model keeps its model methods"
  ),
  stringsAsFactors = FALSE
)

## What a package has under names like `name`: the closest names of the same
## kind, to point at a rename
similar_names <- function(package, name, kind) {
  candidates <- switch(kind,
    `function` = ,
    environment = ls(asNamespace(package), all.names = TRUE),
    private_field = names(cmdstanr:::CmdStanModel$private_fields),
    file = list.files(system.file(dirname(name), package = package))
  )
  target <- basename(name)
  distance <- utils::adist(target, candidates, ignore.case = TRUE)[1, ]
  close <- candidates[order(distance)][seq_len(min(3, length(candidates)))]
  close[distance[order(distance)][seq_along(close)] <= max(3, nchar(target) %/% 3)]
}

## Is one declared internal there?
internal_is_there <- function(package, name, kind) {
  switch(kind,
    `function` = is.function(asNamespace(package)[[name]]),
    environment = is.environment(asNamespace(package)[[name]]),
    file = nzchar(system.file(name, package = package)),
    private_field = name %in% names(cmdstanr:::CmdStanModel$private_fields)
  )
}

## The internals, each checked: character(0) when all are there, else what
## is gone, with the names that look like it (a rename?)
shared_build_problems <- function() {
  problems <- character(0)
  for (i in seq_len(nrow(shared_build_internals))) {
    row <- shared_build_internals[i, ]
    if (internal_is_there(row$package, row$name, row$kind)) next
    similar <- similar_names(row$package, row$name, row$kind)
    problems <- c(problems, sprintf(
      "%s %s %s is gone (%s)%s", row$package, row$kind, row$name, row$use,
      if (length(similar) == 0) "" else sprintf("; renamed? it has %s", toString(similar))
    ))
  }
  # the key, from flatcor's own code: compile_stan_model() names the Bayes
  # factor model paste(backend, "methods")
  if (internal_is_there("flatcor", "compile_stan_model", "function")) {
    code <- paste(deparse(asNamespace("flatcor")$compile_stan_model), collapse = " ")
    key_parts <- c("\"methods\"", "compiled_models[[key]]")
    if (!all(vapply(key_parts, grepl, logical(1), x = code, fixed = TRUE))) {
      problems <- c(problems, sprintf(
        "flatcor's compile_stan_model() no longer keeps the Bayes factor model as \"%s\"",
        shared_build_key
      ))
    }
  }
  problems
}


## ---- Build it, or load it -------------------------------------------------------------

## flatcor's Bayes factor model from the shared build: the executable and
## the model methods, built the first time, loaded every time after
shared_flatcor_model <- function() {
  folder <- shared_build_folder()
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  stan_file <- system.file("stan", "flatcor.stan", package = "flatcor", mustWork = TRUE)

  # the executable: cmdstanr reuses it when it is there and up to date
  model <- cmdstanr::cmdstan_model(stan_file, dir = folder, quiet = TRUE)

  # the model methods, as cmdstanr builds them (expose_model_methods()), but
  # from a file that stays, so that Rcpp finds its build of it next time
  methods_file <- file.path(folder, "flatcor_model_methods.cpp")
  if (!file.exists(methods_file)) {
    hpp <- cmdstanr:::get_standalone_hpp(stan_file, character(0))
    methods <- readLines(system.file("include", "model_methods.cpp", package = "cmdstanr"))
    writeLines(c(hpp, methods), methods_file)
  }
  flags <- function(name) cmdstanr:::get_cmdstan_flags(name)
  libraries <- c("LDLIBS", "LIBSUNDIALS", "TBB_TARGETS", "LDFLAGS_TBB", "SUNDIALS_TARGETS")
  includes <- paste0(" -I\"", system.file("include", package = "cmdstanr"), "\"")
  methods_env <- new.env()
  methods_env$hpp_code_ <- methods_file # cmdstanr only checks that it is set
  tbb <- file.path(cmdstanr::cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb")
  withr::with_path(
    cmdstanr:::repair_path(tbb),
    withr::with_makevars(
      c(
        USE_CXX14 = 1,
        PKG_CPPFLAGS = flags("CPPFLAGS"),
        PKG_CXXFLAGS = paste0(flags("CXXFLAGS"), includes),
        PKG_LIBS = paste(vapply(libraries, flags, ""), collapse = " ")
      ),
      Rcpp::sourceCpp(file = methods_file, env = methods_env, cacheDir = file.path(folder, "rcpp"))
    )
  )
  model$.__enclos_env__$private$model_methods_env_ <- methods_env
  model
}

## Gives flatcor the shared build, so this process's Bayes factors use it;
## when an internal it relies on has changed, flatcor builds its own (the
## 45 seconds), with a warning
use_shared_flatcor_model <- function() {
  problems <- shared_build_problems()
  if (length(problems) > 0) {
    warning(
      "The shared Stan build is not possible (", paste(problems, collapse = "; "),
      "): flatcor builds its own Bayes factor model in this process. ",
      "Test: Rscript R/shared/shared_stan_build.R",
      call. = FALSE
    )
    compile_stan_model <- utils::getFromNamespace("compile_stan_model", "flatcor")
    return(invisible(compile_stan_model("cmdstanr", model_methods = TRUE)))
  }
  model <- shared_flatcor_model()
  assign(shared_build_key, model, envir = asNamespace("flatcor")$compiled_models)
  invisible(model)
}


## ---- The explicit test --------------------------------------------------------------

## One Bayes factor in a fresh R process, with flatcor's own build or the
## shared one: its value, its error and how long it took
bayes_factor_in_fresh_process <- function(shared) {
  callr::r(
    function(shared, project) {
      setwd(project)
      if (!nzchar(Sys.getenv("R_USER_CACHE_DIR"))) {
        Sys.setenv(R_USER_CACHE_DIR = file.path(project, "cache", "r_user_cache"))
      }
      started <- Sys.time()
      if (shared) {
        source(file.path("R", "shared", "shared_stan_build.R"))
        use_shared_flatcor_model()
      }
      set.seed(1)
      x1 <- stats::rnorm(60)
      x2 <- stats::rnorm(60)
      data <- data.frame(y = 0.4 * x1 + stats::rnorm(60), x1, x2)
      result <- suppressWarnings(flatcor::bayes_factor_models(
        y ~ x1 + x2, y ~ x2, data,
        prior = "flat_partial_r",
        iter_sampling = 1000,
        seed = 7
      ))
      c(
        log_bf = result$log_bayes_factor,
        error = result$error,
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    },
    args = list(shared = shared, project = getwd())
  )
}

check_shared_build <- function() {
  problems <- shared_build_problems()
  if (length(problems) > 0) {
    stop("The shared Stan build is not possible:\n  - ", paste(problems, collapse = "\n  - "),
      call. = FALSE
    )
  }
  cat("1. the internals the shared build relies on are there\n")

  own <- bayes_factor_in_fresh_process(shared = FALSE)
  cat(sprintf(
    "2. flatcor's own build:    log BF %.6f, error %.6f, %.0f s\n",
    own[["log_bf"]], own[["error"]], own[["seconds"]]
  ))
  shared_first <- bayes_factor_in_fresh_process(shared = TRUE)
  shared_again <- bayes_factor_in_fresh_process(shared = TRUE)
  cat(sprintf(
    "3. the shared build:       log BF %.6f, error %.6f, %.0f s (then %.0f s)\n",
    shared_again[["log_bf"]], shared_again[["error"]],
    shared_first[["seconds"]], shared_again[["seconds"]]
  ))

  values <- function(run) unname(run[c("log_bf", "error")])
  same <- identical(values(own), values(shared_first)) &&
    identical(values(own), values(shared_again))
  if (!same) {
    stop("The shared build gives another Bayes factor than flatcor's own build", call. = FALSE)
  }
  cat("The shared build gives the same Bayes factor as flatcor's own build.\n")
  invisible(TRUE)
}

## Run as a script (Rscript R/shared/shared_stan_build.R): the explicit test
if (sys.nframe() == 0) check_shared_build()
