## =====================================================================
##  Does BridgeStan give flatcor's Bayes factors, on Linux in the image?
##
##  flatcor's Bayes factors need the Stan model's log density in R. Now that
##  comes from cmdstanr's model methods, built again in every R process (about
##  45 s per worker per run). BridgeStan (from the Stan developers) compiles
##  the log density into a library file any process loads. On a Mac the two
##  gave bit-identical Bayes factors; this checks the same on Linux, inside
##  the validation's image (.github/workflows/bridgestan.yml):
##
##    Rscript ci/bridgestan_check.R
##
##    1. a Bayes factor from flatcor as it is and from the BridgeStan route
##       (the same data steps, sampler settings, seeds and bridgesampling
##       call; only the log density from BridgeStan), for three model pairs
##    2. the same Bayes factor in four fresh R processes at once
##    3. 50 Bayes factors in a row in one process (a long-lived worker):
##       time and memory
##
##  Nothing here changes flatcor or the pipeline. It stops with an error when a
##  result is not a number, when the two routes differ by more than 1e-8 on
##  the log scale (far below a Bayes factor's Monte Carlo error, about 0.01),
##  or when the four processes disagree; whether the routes are identical is
##  reported either way.
## =====================================================================

## ---- BridgeStan: the R package and its C++ sources, in a temporary folder ------------
## (not in the image: this is a check, not yet part of the environment)
setup <- file.path(tempdir(), "bridgestan")
dir.create(file.path(setup, "lib"), recursive = TRUE, showWarnings = FALSE)
install.packages(
  "bridgestan",
  lib = file.path(setup, "lib"),
  repos = c("https://r-multiverse.r-universe.dev", "https://cloud.r-project.org"),
  quiet = TRUE
)
library(bridgestan, lib.loc = file.path(setup, "lib"))
version <- as.character(utils::packageVersion("bridgestan", lib.loc = file.path(setup, "lib")))
sources <- file.path(setup, sprintf("bridgestan-%s.tar.gz", version))
utils::download.file(
  sprintf(
    "https://github.com/roualdes/bridgestan/releases/download/v%s/bridgestan-%s.tar.gz",
    version, version
  ),
  sources,
  quiet = TRUE
)
utils::untar(sources, exdir = setup)
bridgestan_path <- file.path(setup, sprintf("bridgestan-%s", version))
Sys.setenv(BRIDGESTAN = bridgestan_path)

## flatcor's Stan program, compiled by BridgeStan into a library file
stan_copy <- file.path(setup, "model", "flatcor.stan")
dir.create(dirname(stan_copy), showWarnings = FALSE)
file.copy(system.file("stan", "flatcor.stan", package = "flatcor"), stan_copy)
library_file <- sub("[.]stan$", "_model.so", stan_copy)
started <- Sys.time()
# make's tools (ar) write temporary files to TMPDIR; when that names a folder
# that does not exist (R itself then falls back quietly), give them R's own
temporary <- Sys.getenv("TMPDIR")
if (!nzchar(temporary) || !dir.exists(temporary)) temporary <- tempdir()
status <- withr::with_envvar(
  c(TMPDIR = temporary),
  system2("make", c("-C", bridgestan_path, library_file))
)
if (status != 0 || !file.exists(library_file)) stop("BridgeStan could not compile flatcor.stan")
cat(sprintf(
  "BridgeStan %s compiled flatcor.stan in %.0f s (with BridgeStan's own first build)\n",
  version, as.numeric(difftime(Sys.time(), started, units = "secs"))
))


## ---- One Bayes factor by the BridgeStan route ---------------------------------------
## As flatcor::bayes_factor_models(): its data steps, sampler settings, seeds and
## bridgesampling call, with fitting by cmdstanr's plain build and the log
## density from BridgeStan (the same density as cmdstanr's log_prob, which
## flatcor uses now: up to a constant, with the Jacobian)
bridgestan_bayes_factor <- function(formulas, data, prior, seed, library_file,
                                    iter_sampling = 1000) {
  plain <- utils::getFromNamespace("compile_stan_model", "flatcor")("cmdstanr")
  rows <- utils::getFromNamespace("complete_rows", "flatcor")(data, formulas)
  per_model <- lapply(c(a = 0, b = 1), function(offset) {
    parts <- utils::getFromNamespace("predictors_and_outcome", "flatcor")(
      formulas[[offset + 1]], rows
    )
    k <- ncol(parts$predictors)
    stan_data <- list(
      N = length(parts$outcome), K = k, X = parts$predictors, y = parts$outcome,
      prior_only = 0L, proper_stand_ins = 0L, per_observation = 0L,
      r_squared_shape = utils::getFromNamespace("prior_shape", "flatcor")(prior, k)
    )
    fit <- plain$sample(
      data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000,
      iter_sampling = iter_sampling, adapt_delta = 0.95, seed = seed + offset,
      refresh = 0, show_messages = FALSE
    )

    # the log density: BridgeStan's library, loaded with this model's data
    data_file <- tempfile(fileext = ".json")
    cmdstanr::write_stan_json(stan_data, data_file)
    density <- bridgestan::StanModel$new(library_file, data_file, seed = 1, warn = FALSE)
    unlink(data_file)

    # the draws on Stan's unconstrained scale (BridgeStan names elements z.1,
    # cmdstanr z[1])
    names_in_fit <- sub("^([^.]+)[.]([0-9]+)$", "\\1[\\2]", density$param_names())
    draws <- posterior::as_draws_matrix(fit$draws(names_in_fit))
    unconstrained <- t(apply(unclass(draws), 1, density$param_unconstrain))
    colnames(unconstrained) <- density$param_unc_names()

    # the argument names bridgesampling's own log_posterior functions use
    log_posterior <- function(s.row, data) { # nolint: object_name_linter.
      tryCatch(
        data$log_density(as.numeric(s.row), propto = TRUE, jacobian = TRUE),
        error = function(e) -Inf
      )
    }
    bound <- stats::setNames(rep(Inf, ncol(unconstrained)), colnames(unconstrained))
    with_seed <- utils::getFromNamespace("with_seed", "flatcor")
    bridge <- with_seed(seed + offset, bridgesampling::bridge_sampler(
      unconstrained,
      log_posterior = log_posterior,
      data = density,
      lb = -bound,
      ub = bound,
      method = "normal",
      silent = TRUE
    ))
    c(logml = bridge$logml, cv = bridgesampling::error_measures(bridge)$cv)
  })
  c(
    log_bf = per_model$a[["logml"]] - per_model$b[["logml"]],
    error = sqrt(per_model$a[["cv"]]^2 + per_model$b[["cv"]]^2)
  )
}

## A dataset of n rows: x1 with a slope, x2 and x3 without
dataset <- function(data_seed, n = 60) {
  set.seed(data_seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rnorm(n)
  data.frame(y = 0.4 * x1 + stats::rnorm(n), x1, x2, x3)
}

problems <- character(0)


## ---- 1. The two routes, three model pairs ---------------------------------------------
pairs <- list(
  list(formulas = list(y ~ x1 + x2, y ~ x2), prior = "flat_partial_r"),
  list(formulas = list(y ~ x1 + x2, y ~ 1), prior = "flat_r"),
  list(formulas = list(y ~ x1 + x2 + x3, y ~ x2 + x3), prior = "flat_r_squared")
)
cat("\n1. flatcor as it is against the BridgeStan route\n")
for (pair in pairs) {
  data <- dataset(1)
  label <- sprintf(
    "%s vs %s, %s", deparse(pair$formulas[[1]][[3]]), deparse(pair$formulas[[2]][[3]]), pair$prior
  )
  own <- suppressWarnings(flatcor::bayes_factor_models(
    pair$formulas[[1]], pair$formulas[[2]], data,
    prior = pair$prior,
    iter_sampling = 1000,
    seed = 7
  ))
  route <- bridgestan_bayes_factor(pair$formulas, data, pair$prior, seed = 7, library_file)
  difference <- route[["log_bf"]] - own$log_bayes_factor
  cat(sprintf(
    paste0(
      "   %s\n",
      "     flatcor:    log BF %.17g, error %.17g\n",
      "     BridgeStan: log BF %.17g, error %.17g\n",
      "     identical: %s (log BF), %s (error); difference %.3g\n"
    ),
    label, own$log_bayes_factor, own$error, route[["log_bf"]], route[["error"]],
    identical(route[["log_bf"]], own$log_bayes_factor), identical(route[["error"]], own$error),
    difference
  ))
  if (!all(is.finite(c(route, own$log_bayes_factor)))) {
    problems <- c(problems, paste(label, ": not a number"))
  } else if (abs(difference) > 1e-8) {
    problems <- c(problems, sprintf("%s: the routes differ by %.3g", label, difference))
  }
}


## ---- 2. Four fresh R processes at once ----------------------------------------------
cat("\n2. Four fresh R processes at once, each loading the library\n")
check_script <- normalizePath("ci/bridgestan_check.R")
processes <- lapply(1:4, function(i) {
  callr::r_bg(
    function(script, library_file, bridgestan_lib) {
      library(bridgestan, lib.loc = bridgestan_lib)
      # the functions of this script, without running it
      code <- parse(script)
      keep <- vapply(code, function(e) {
        is.call(e) && identical(e[[1]], as.name("<-")) &&
          as.character(e[[2]]) %in% c("bridgestan_bayes_factor", "dataset")
      }, logical(1))
      for (e in code[keep]) eval(e, globalenv())
      started <- Sys.time()
      result <- bridgestan_bayes_factor(
        list(y ~ x1 + x2, y ~ x2), dataset(1), "flat_partial_r", seed = 7, library_file
      )
      c(result, seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    },
    args = list(
      script = check_script, library_file = library_file,
      bridgestan_lib = file.path(setup, "lib")
    )
  )
})
for (p in processes) p$wait()
fresh <- lapply(processes, function(p) p$get_result())
for (i in seq_along(fresh)) {
  cat(sprintf(
    "   process %d: log BF %.17g, %.1f s\n", i, fresh[[i]][["log_bf"]], fresh[[i]][["seconds"]]
  ))
}
if (length(unique(vapply(fresh, `[[`, numeric(1), "log_bf"))) != 1) {
  problems <- c(problems, "the four processes disagree")
}


## ---- 3. 50 Bayes factors in a row in one process --------------------------------------
cat("\n3. 50 Bayes factors in a row in one process\n")
## this process's memory (RSS): from /proc on Linux, from ps elsewhere
memory_mb <- function() {
  if (file.exists("/proc/self/status")) {
    status <- readLines("/proc/self/status")
    return(as.numeric(gsub("[^0-9]", "", grep("^VmRSS", status, value = TRUE))) / 1024)
  }
  as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE)) / 1024
}
priors <- c("flat_partial_r", "flat_r", "flat_r_squared")
in_a_row <- do.call(rbind, lapply(1:50, function(i) {
  started <- Sys.time()
  result <- bridgestan_bayes_factor(
    list(y ~ x1 + x2, y ~ x2), dataset(100 + i), priors[(i - 1) %% 3 + 1],
    seed = 1000 + i, library_file
  )
  gc()
  data.frame(
    log_bf = result[["log_bf"]], error = result[["error"]],
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    memory_mb = memory_mb()
  )
}))
cat(sprintf(
  paste0(
    "   all numbers: %s; seconds per Bayes factor: first 10 %.2f, last 10 %.2f;",
    " memory (MB) after 1, 10, 50: %.0f, %.0f, %.0f\n"
  ),
  all(is.finite(in_a_row$log_bf) & is.finite(in_a_row$error)),
  mean(in_a_row$seconds[1:10]), mean(in_a_row$seconds[41:50]),
  in_a_row$memory_mb[1], in_a_row$memory_mb[10], in_a_row$memory_mb[50]
))
if (!all(is.finite(in_a_row$log_bf) & is.finite(in_a_row$error))) {
  problems <- c(problems, "a Bayes factor in the row of 50 is not a number")
}


## ---- The verdict ----------------------------------------------------------------------
if (length(problems) > 0) {
  stop("BridgeStan check:\n  - ", paste(problems, collapse = "\n  - "), call. = FALSE)
}
cat("\nBridgeStan gives flatcor's Bayes factors here.\n")
