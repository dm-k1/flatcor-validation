## =====================================================================
##  The reference task: one fixed flatcor Bayes factor, timed
##
##    Rscript R/timing/reference_task.R
##
##  prints how many seconds this machine takes for it: the median of three
##  runs, after one untimed run that compiles the model. Every timing run
##  (time_study.R) records it on its runner, so the estimate
##  (estimate_run_time.R) can count each timing in multiples of it: a faster
##  or slower processor, or a busy host, then divides out. Run it on the
##  machine you plan to use and give its seconds to the estimate
##  (--reference=<seconds>) to read the estimate off for that machine.
##
##  The task: n = 1000, three predictors against two, flat partial r, the
##  default draws (4 chains of 10,000, run in parallel), fixed data and seed.
##  Its seeds stay fixed on purpose and do not follow root_seed
##  (R/shared/random_streams.R): it is a yardstick for machines, so it must
##  be the same task on every machine and in every run.
## =====================================================================

reference_data <- function() {
  set.seed(2026)
  data <- data.frame(x1 = rnorm(1000), x2 = rnorm(1000), x3 = rnorm(1000))
  data$y <- 0.3 * data$x1 + 0.2 * data$x2 + rnorm(1000)
  data
}

reference_run <- function(data) {
  started <- Sys.time()
  suppressMessages(flatcor::bayes_factor_models(
    y ~ x1 + x2 + x3,
    y ~ x2 + x3,
    data,
    prior = "flat_partial_r",
    seed = 1
  ))
  as.numeric(difftime(Sys.time(), started, units = "secs"))
}

reference_seconds <- function() {
  data <- reference_data()
  reference_run(data) # untimed: compiles the model
  stats::median(vapply(1:3, function(i) reference_run(data), numeric(1)))
}

## run as a script (not sourced by time_study.R): print the seconds
if (sys.nframe() == 0) {
  cat(sprintf(
    "Reference task: %.1f seconds on this machine (%d cores)\n",
    reference_seconds(), parallelly::availableCores()
  ))
}
