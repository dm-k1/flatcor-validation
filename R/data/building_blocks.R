## =====================================================================
##  Building blocks of every data set (docs/validation_tests.md, "The data
##  sets")
##
##  Every dataset has an ID: its data set, its cell and its replication. Its
##  random numbers come from its own random stream, seeded from that ID, so
##  the same ID always gives the same data, in any study, on any machine, in
##  any order, on any number of cores. Datasets are made again from their
##  ID, never stored.
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

## n x k standard normal predictors with the same correlation between every
## pair, named x1 ... xk; centred if asked (the model conditions on them)
correlated_predictors <- function(n, k, correlation, centred = FALSE) {
  population <- matrix(correlation, k, k)
  diag(population) <- 1

  predictors <- matrix(rnorm(n * k), n) %*% chol(population)
  if (centred) predictors <- sweep(predictors, 2, colMeans(predictors))
  colnames(predictors) <- paste0("x", seq_len(k))
  predictors
}
