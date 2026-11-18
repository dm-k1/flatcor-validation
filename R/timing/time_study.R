## =====================================================================
##  Times one study on a fraction of its datasets
##
##    VALIDATION_TIMING=0.04 Rscript R/timing/time_study.R 3a
##
##  runs the study's steps of the pipeline (_targets.R: the targets ending in
##  _3a) in timing mode (R/shared/pipeline.R): the full design, but 4% of the
##  datasets (at least 2), in its own store under cache/timing/, so it never
##  mixes with a real run. It writes cache/timing/time_3a_at_0.04.rds with the
##  elapsed time and the cores. R/timing/estimate_run_time.R scales these up to
##  the full run; the "Run-time estimate" workflow (.github/workflows/timing.yml)
##  runs both. Studies 1a (it does not depend on the number of datasets) and
##  1b (20 datasets per prior, too few to scale) are timed at fraction 1.
## =====================================================================
study <- commandArgs(trailingOnly = TRUE)[1]
fraction <- as.numeric(Sys.getenv("VALIDATION_TIMING", NA))
stopifnot(
  "give a study, as in Rscript R/timing/time_study.R 3a" = !is.na(study),
  "set VALIDATION_TIMING to the fraction of datasets, as in 0.04" = !is.na(fraction)
)

## the reference task first, on this machine (R/timing/reference_task.R), so
## the estimate can divide out how fast this machine is
source(file.path("R", "timing", "reference_task.R"))
reference <- tryCatch(reference_seconds(), error = function(e) NA_real_)

store <- file.path("cache", "timing", paste0("store_", study, "_", fraction))
started <- Sys.time()
error <- tryCatch(
  {
    targets::tar_make(
      names = tidyselect::ends_with(!!paste0("_", study)),
      store = store,
      reporter = "silent"
    )
    NA_character_
  },
  error = function(e) conditionMessage(e)
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

dir.create(file.path("cache", "timing"), recursive = TRUE, showWarnings = FALSE)
cpu_model <- tryCatch(
  {
    model_line <- grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
    sub(".*:\\s*", "", model_line)
  },
  error = function(e) NA_character_,
  warning = function(w) NA_character_
)
saveRDS(
  list(
    study = study,
    fraction = fraction,
    cores = parallelly::availableCores(),
    cpu_model = cpu_model,
    elapsed_seconds = elapsed,
    reference_seconds = reference,
    error = error
  ),
  file.path("cache", "timing", sprintf("time_%s_at_%s.rds", study, fraction))
)
cat(sprintf(
  "\nStudy %s at fraction %s: %.1f minutes elapsed%s\n", study, fraction, elapsed / 60,
  if (is.na(error)) "" else paste0("; stopped with an error: ", error)
))
if (!is.na(error)) quit(status = 1)
