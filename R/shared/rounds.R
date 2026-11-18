## =====================================================================
##  Rounds: how many datasets each cell gets
##  (docs/validation_tests.md, "How many datasets")
##
##  Every quantity has an accuracy Δ; its Monte Carlo SE must be at most
##  Δ / 1.96. A study that measures its numbers runs in rounds (the
##  Chow-Robbins rule): every cell starts with a first round, fixed in
##  advance; after each round, a cell whose quantities have not all reached
##  their target gets the datasets it is estimated to lack, and the round
##  is run again. The plan, how many datasets each cell has so far, is a
##  small file in the targets store, so a fresh start deletes it with the
##  store. run_all.R runs the rounds; the pipeline only reads the plan.
## =====================================================================


## ---- The target ------------------------------------------------------------------------

## The largest Monte Carlo SE at which a deviation of delta reaches the 95% line
target_mcse <- function(delta) delta / 1.96

## Δ, the target and whether it is reached, added to a results table with
## columns `quantity` and `mcse`; `deltas`: Δ per quantity, named as in the
## table, or NULL when the table has a column `delta` already (a Δ that
## differs per row, such as 10% of each row's expected value)
with_targets <- function(table, deltas = NULL) {
  if (!is.null(deltas)) table$delta <- unname(deltas[table$quantity])
  table$target_mcse <- target_mcse(table$delta)
  # a quantity without a Δ (an overall check) has no target to wait for
  table$reached <- is.na(table$delta) | table$mcse <= table$target_mcse
  table
}


## ---- The plan ------------------------------------------------------------------------

## An .rds file, not .csv: it keeps every column's type, and a piece of work
## is the same piece only if its cell's values are identical (2 is not 2L)
plan_file <- function(study, store) file.path(store, "rounds", paste0(study, ".rds"))

## The plan of a study: its cells, each with the number of datasets so far.
## Without a plan file (the first round, or a timing run), every cell has
## the first round.
read_plan <- function(study, store, cells, first_round) {
  file <- plan_file(study, store)
  if (file.exists(file)) {
    return(readRDS(file))
  }
  data.frame(cells, datasets = first_round, stringsAsFactors = FALSE)
}

write_plan <- function(plan, study, store) {
  dir.create(dirname(plan_file(study, store)), recursive = TRUE, showWarnings = FALSE)
  saveRDS(plan, plan_file(study, store))
}

## The pieces of work of a plan: one row per cell and batch, with the cell's
## columns and the first and last replication of the batch. Batches are cut
## from 1 in steps of batch_size, so a later round only adds rows: the
## pieces already done stay done.
pieces_of_plan <- function(plan, batch_size) {
  rows <- lapply(seq_len(nrow(plan)), function(i) {
    starts <- seq(1, plan$datasets[i], by = batch_size)
    cell <- plan[rep(i, length(starts)), setdiff(names(plan), "datasets"), drop = FALSE]
    cell$from <- starts
    cell$to <- pmin(starts + batch_size - 1, plan$datasets[i])
    cell
  })
  pieces <- do.call(rbind, rows)
  rownames(pieces) <- NULL
  pieces
}


## ---- The next round ----------------------------------------------------------------------

## The plan for the next round: every cell with a quantity whose Monte Carlo
## SE is above its target gets the datasets it is estimated to lack (the SE
## shrinks with the square root of the datasets: it needs
## datasets x (mcse / target)^2), rounded up to whole batches, at least one.
## `table`: the study's results table with_targets(); `by`: the columns that
## name the unit a quantity is computed over, in the plan and in the table:
## a cell, or several (a slope over n: then every cell of the unit, each n,
## gets the datasets). `one_batch`: add one batch
## only (the smoke test, which tests the step, not the numbers). Returns the
## plan unchanged when every cell has reached its targets.
next_plan <- function(plan, table, by, batch_size, one_batch = FALSE) {
  shortfall <- table[!table$reached, c(by, "mcse", "target_mcse"), drop = FALSE]
  if (nrow(shortfall) == 0) {
    return(plan)
  }
  shortfall$needed_share <- (shortfall$mcse / shortfall$target_mcse)^2
  worst <- stats::aggregate(needed_share ~ ., data = shortfall[c(by, "needed_share")], FUN = max)

  plan_key <- do.call(paste, plan[by])
  for (i in seq_len(nrow(worst))) {
    rows <- which(plan_key == do.call(paste, worst[i, by, drop = FALSE]))
    lacking <- plan$datasets[rows] * (worst$needed_share[i] - 1)
    batches <- if (one_batch) 1 else pmax(1, ceiling(lacking / batch_size))
    plan$datasets[rows] <- plan$datasets[rows] + batch_size * batches
  }
  plan
}
