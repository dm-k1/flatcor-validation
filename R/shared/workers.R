## =====================================================================
##  What every worker starts with: its environment and the validation's
##  code. Sourced by _targets.R (the pipeline's workers) and run_all.R (the
##  queue's workers, R/shared/queue.R) before any worker starts, so both
##  kinds inherit the same settings.
## =====================================================================

## flatcor compiles its Stan model into R's user cache folder: keep that in
## cache/, writable everywhere
if (!nzchar(Sys.getenv("R_USER_CACHE_DIR"))) {
  Sys.setenv(R_USER_CACHE_DIR = file.path(getwd(), "cache", "r_user_cache"))
}

## Linux: every Stan model loads CmdStan's TBB library (libtbb.so.2) when it
## starts, and compiling a model can make CmdStan rebuild that library in
## place. With several workers, one worker's compile then removes it under a
## model another worker is starting ("libtbb.so.2: cannot open shared object
## file", seen with 7 workers on AWS). A copy made here, before any worker
## starts, and put first on the library path the workers inherit, is never
## touched by such a rebuild.
if (Sys.info()[["sysname"]] == "Linux") {
  tbb_copy <- file.path(getwd(), "cache", "tbb")
  dir.create(tbb_copy, recursive = TRUE, showWarnings = FALSE)
  tbb_libraries <- list.files(
    file.path(cmdstanr::cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb"),
    pattern = "^libtbb.*[.]so[.]2$", full.names = TRUE
  )
  file.copy(tbb_libraries, tbb_copy, overwrite = TRUE)
  library_path <- Sys.getenv("LD_LIBRARY_PATH")
  if (!startsWith(library_path, tbb_copy)) {
    library_path <- paste(c(tbb_copy, library_path[nzchar(library_path)]), collapse = ":")
    Sys.setenv(LD_LIBRARY_PATH = library_path)
  }
}

## One thread per process for linear algebra (BLAS, OpenMP). Each piece of
## work already keeps 4 cores busy with Stan's chains, and the workers are
## sized for that; a BLAS that starts a thread per core in every worker
## overloads the machine instead (a load of about 200 on 64 cores, seen on
## AWS). A value set from outside (a cluster job, say) is kept.
local({
  for (variable in c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")) {
    if (!nzchar(Sys.getenv(variable))) do.call(Sys.setenv, stats::setNames(list("1"), variable))
  }
})

## Each piece of work runs Stan's 4 chains in parallel, so a quarter of the
## cores work on pieces at the same time, each in its own R process. The
## cores are those given to this run: on a cluster, the job's allocation
## (parallelly reads SLURM, PBS and cgroups), not all cores of the node.
worker_count <- function() max(1, (parallelly::availableCores() - 1) %/% 4)

## The validation's code, in the order it is sourced: the shared code, then
## each study's design and code, and last, only to record or replay the
## expensive calls, the dry run's file, which replaces a few of the
## functions before it (R/dry_run/replay.R)
code_files <- function() {
  studies <- c(
    "1a_prior_sampling", "1b_invariance", "2a_posterior",
    "2b_estimation", "3a_bayes_factors", "3b_default_criteria"
  )
  replay <- identical(Sys.getenv("VALIDATION_DRY"), "1") ||
    identical(Sys.getenv("VALIDATION_RECORD"), "1")
  c(
    file.path("R", "shared", c(
      "random_streams.R", "pipeline.R", "expectations.R", "results.R", "per_dataset.R",
      "comparison_methods.R", "rounds.R"
    )),
    list.files(file.path("R", "1_data"), pattern = "[.]R$", full.names = TRUE),
    file.path("R", "2_studies", studies, "design.R"),
    file.path("R", "2_studies", studies, "code.R"),
    if (replay) file.path("R", "dry_run", "replay.R")
  )
}
