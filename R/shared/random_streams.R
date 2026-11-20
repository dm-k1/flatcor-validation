## =====================================================================
##  Randomness: every random number of the validation follows from one
##  number, root_seed, and a name (docs/validation_tests.md, "Randomness").
##
##  Everything random has a name that says what it is: a dataset's name is
##  its data set, cell and replication; a Stan fit's name is its dataset's
##  name and the prior; a bootstrap's name is what it bootstraps. The name,
##  hashed with root_seed, gives a seed:
##    with_stream_of(name, code)  seeds R's own random numbers (data,
##                                bootstraps, BayesFactor's draws)
##    seed_of(...)                a seed for Stan and bridge sampling
##  So the same name always gives the same numbers, in any study, on any
##  machine, in any order, on any number of cores, in any round. Datasets
##  are made again from their name, never stored.
## =====================================================================


## The one seed of the validation. Change it (and commit) to run the whole
## validation again with fresh random numbers: every step then reruns.
## targets' own seed follows it too (_targets.R), and every results file
## records it. With 0, every number is the same as in v1.
root_seed <- 0L


## Evaluate `code` with the random stream of one name, and restore the
## caller's random state afterwards, so making a dataset never changes what
## the rest of a study draws. The seed is a hash of the name, so
## neighbouring names do not get neighbouring seeds.
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
  set.seed(digest::digest2int(dataset_id, seed = root_seed))
  code
}

## A seed for Stan and bridge sampling, from a name (a dataset's ID and a
## label): the same hash as with_stream_of(), as a positive whole number.
## flatcor gives the second model of a Bayes factor seed + 1.
seed_of <- function(...) {
  abs(digest::digest2int(paste(...), seed = root_seed)) %% .Machine$integer.max
}
