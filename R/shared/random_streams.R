## =====================================================================
##  Random streams: every dataset, and every other random draw of a study,
##  has an ID (for a dataset: its data set, its cell and its replication).
##  Its random numbers come from its own random stream, seeded from that ID,
##  so the same ID always gives the same numbers, in any study, on any
##  machine, in any order, on any number of cores. Datasets are made again
##  from their ID, never stored (docs/validation_tests.md, "The data sets").
##
##  In short: name -> hash -> seed. This file seeds R's own random numbers
##  (data, bootstraps, BayesFactor's draws); seed_of() in
##  comparison_methods.R makes the seed for Stan and bridge sampling from a
##  name in the same way. The hash's own seed is left at 0: it is the one
##  number from which every random number of the validation follows.
## =====================================================================


## Evaluate `code` with the random stream of one dataset ID, and restore the
## caller's random state afterwards, so making a dataset never changes what
## the rest of a study draws. The seed is a hash of the ID, so neighbouring
## IDs do not get neighbouring seeds.
with_stream_of <- function(dataset_id, code) {
  caller_kind <- RNGkind()
  caller_seed <- get0(".Random.seed", envir = globalenv())
  on.exit({
    do.call(RNGkind, as.list(caller_kind))
    if (is.null(caller_seed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", caller_seed, envir = globalenv())
    }
  })

  RNGkind("L'Ecuyer-CMRG")
  set.seed(digest::digest2int(dataset_id))
  code
}
