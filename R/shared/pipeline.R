## =====================================================================
##  Helpers the studies share in the pipeline (_targets.R)
## =====================================================================

## ---- How much to run ---------------------------------------------------------------
## Full run by default. VALIDATION_SMOKE=1: a corner of each design, to test
## the code (results go to cache/, never mixed with real ones).
## VALIDATION_TIMING=<fraction>: the full design with that fraction of the
## datasets (at least 2), to estimate the run time. The floor is 2, not
## more, so that the two fractions of the estimate (4% and 6%) give different
## numbers of datasets also for study 3b (50 or 100 per cell: 2 and 3, or 4
## and 6), and both stay small enough for a 6-hour GitHub runner.
smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
timing_fraction <- suppressWarnings(as.numeric(Sys.getenv("VALIDATION_TIMING", NA)))

## The number of datasets to run: `full` normally, `in_smoke` in the smoke test
datasets <- function(full, in_smoke) {
  if (smoke) {
    return(in_smoke)
  }
  if (!is.na(timing_fraction)) {
    return(min(full, max(2, ceiling(full * timing_fraction))))
  }
  full
}

## The replications 1 ... total in batches of `size`: each batch is one piece
## of work, saved as soon as it is done
split_into_batches <- function(total, size) split(seq_len(total), ceiling(seq_len(total) / size))

## ---- flatcor's Stan models, compiled one worker at a time ---------------------------
## Compiling a model can make CmdStan rebuild shared files in its folder (the
## main.o every model is linked with, the TBB library), so two compiles at
## once break each other ("cannot find src/cmdstan/main.o", seen with 7
## workers on AWS). Each worker compiles flatcor's two versions of the model
## (for estimation, and with model methods for Bayes factors) once, here,
## under a lock: one worker at a time, the others wait. flatcor keeps them in
## memory, so the worker's later fits and Bayes factors compile nothing. The
## lock is a folder: creating one is atomic, so only one worker succeeds; a
## lock older than 30 minutes is taken to be left by a worker that stopped
## while compiling, and removed.
compile_flatcor_models <- function() {
  if (isTRUE(getOption("validation_flatcor_compiled"))) {
    return(invisible(TRUE))
  }
  dir.create("cache", showWarnings = FALSE)
  lock <- file.path("cache", "compile.lock")
  repeat {
    if (dir.create(lock, showWarnings = FALSE)) break
    age <- difftime(Sys.time(), file.mtime(lock), units = "mins")
    if (!is.na(age) && age > 30) unlink(lock, recursive = TRUE)
    Sys.sleep(5)
  }
  on.exit(unlink(lock, recursive = TRUE))

  compile_stan_model <- utils::getFromNamespace("compile_stan_model", "flatcor")
  compile_stan_model("cmdstanr")
  compile_stan_model("cmdstanr", model_methods = TRUE)
  options(validation_flatcor_compiled = TRUE)
  invisible(TRUE)
}

## The model for the studies that sample it directly (the prior only, SBC):
## the pipeline's first step, so that their workers find it compiled
flatcor_stan_model <- function() {
  compile_flatcor_models()
  compile_stan_model <- utils::getFromNamespace("compile_stan_model", "flatcor")
  compile_stan_model("cmdstanr")
}
prior_shape <- function(prior, k) utils::getFromNamespace("prior_shape", "flatcor")(prior, k)

## ---- Monte Carlo SEs ---------------------------------------------------------------------
## Of an estimated proportion. It uses the Agresti-Coull adjusted proportion
## (x + 2) / (n + 4), so a share of 0 or 1 does not get an SE of zero.
mcse_of_proportion <- function(p, n_sim) {
  adjusted <- (p * n_sim + 2) / (n_sim + 4)
  sqrt(adjusted * (1 - adjusted) / (n_sim + 4))
}
