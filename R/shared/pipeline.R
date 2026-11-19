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

## ---- flatcor's Stan model and log density, built once ----------------------------
## flatcor builds two things (flatcor 2.1): the model, for sampling
## (flatcor_model(), cmdstanr's executable) and its log density, for the Bayes
## factors' bridge sampling (flatcor_density(), BridgeStan's library). Both
## are files in flatcor's cache folder, which _targets.R and run_all.R put in
## cache/r_user_cache, so every worker finds them. run_all.R calls this once
## before any worker starts, so the workers only load them. Compiling can make
## CmdStan rebuild files every model shares (the main.o every model is linked
## with, the TBB library), so two compiles at once break each other ("cannot
## find src/cmdstan/main.o", seen with 7 workers on AWS): a worker that still
## has to build does so under a lock, one worker at a time (flatcor also locks
## its own cache folder). The lock is a folder: creating one is atomic, so
## only one worker succeeds; a lock older than 30 minutes is taken to be left
## by a worker that stopped while compiling, and removed.
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

  flatcor::flatcor_model("cmdstanr")
  flatcor::flatcor_density()
  options(validation_flatcor_compiled = TRUE)
  invisible(TRUE)
}

## The model for the studies that sample it directly (the prior only, SBC):
## the pipeline's first step, so that their workers find it compiled
flatcor_stan_model <- function() {
  compile_flatcor_models()
  flatcor::flatcor_model("cmdstanr")
}
prior_shape <- function(prior, k) utils::getFromNamespace("prior_shape", "flatcor")(prior, k)

## ---- Monte Carlo SEs ---------------------------------------------------------------------
## Of an estimated proportion. It uses the Agresti-Coull adjusted proportion
## (x + 2) / (n + 4), so a share of 0 or 1 does not get an SE of zero.
mcse_of_proportion <- function(p, n_sim) {
  adjusted <- (p * n_sim + 2) / (n_sim + 4)
  sqrt(adjusted * (1 - adjusted) / (n_sim + 4))
}
