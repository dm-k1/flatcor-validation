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

## ---- flatcor's Stan model --------------------------------------------------------------
## Compiled once (the pipeline's first step for every study that samples), so
## that parallel workers find the executable instead of compiling it at the
## same time. Studies that sample the model directly (the prior only, SBC) use
## this; fit_regression() and bayes_factor_models() compile it themselves.
flatcor_stan_model <- function() {
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
