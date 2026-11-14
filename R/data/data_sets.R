## =====================================================================
##  The shared data sets of the validation (docs/appendix/C_study_design.md)
##
##  Every dataset has an ID: the data set, its design cell and its
##  replication. Its random numbers come from its own L'Ecuyer-CMRG stream,
##  seeded from that ID, so the same ID always gives the same data, in any
##  study, on any machine, in any order, on any number of cores. Datasets are
##  regenerated from their ID, not stored.
## =====================================================================

## Evaluate code with the random stream of one dataset, given its ID (a
## character string), and restore the caller's random state afterwards, so
## generating a dataset never changes what the rest of a study draws. The seed
## is a hash of the ID, so neighbouring IDs do not get neighbouring seeds.
with_stream_of <- function(dataset_id, code) {
  caller_kind <- RNGkind()
  caller_seed <- get0(".Random.seed", envir = globalenv())
  on.exit({
    do.call(RNGkind, as.list(caller_kind))
    if (is.null(caller_seed)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", caller_seed, envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(digest::digest2int(dataset_id))
  code
}
