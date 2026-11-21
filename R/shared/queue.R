## =====================================================================
##  The queue: the studies that run in rounds, cell by cell
##  (run_all.R runs it; the rule itself is in R/shared/rounds.R)
##
##  Every cell goes through its own rounds. (A "cell" here is what a
##  quantity is computed over: one cell of the design, or the cells of a
##  slope over n together.) Its pieces of work go into the queue; as soon as
##  the last of them is back, the cell is judged on its own results, and if
##  a quantity is short of its target, its next batch goes into the queue
##  at once. No cell waits for another, so the workers stay busy until the
##  last cell is done.
##
##  The results are the same as with rounds over all cells at once: a
##  dataset's numbers come from its name (R/shared/random_streams.R), and a
##  cell's verdict only from its own results, so when a cell is judged
##  changes nothing.
##
##  The budget: a batch that would take the run past it is not started. The
##  cell is then reported as short of its target, with what its next batch
##  would cost; the researcher decides (raise the budget and start again,
##  or report the cell as short). docs/validation_tests.md, "How many datasets".
##
##  Each finished piece is saved in the store (R/shared/rounds.R,
##  run_piece()), so a run that stops continues where it was, and the
##  pipeline reads the pieces afterwards instead of computing them again.
## =====================================================================


## ---- Running pieces ---------------------------------------------------------------

## What a queue worker does with a piece: it loads the validation's code the
## first time (R/shared/workers.R, code_files()) and then runs the piece
## (run_piece(), R/shared/rounds.R). `task`: list(design, piece, store, files).
run_task <- function(task) {
  if (!exists("run_piece", envir = globalenv())) {
    for (file in task$files) sys.source(file, envir = globalenv())
  }
  if (!exists("queue_model", envir = globalenv())) {
    assign("queue_model", flatcor_stan_model(), envir = globalenv())
  }
  run_piece(task$design, task$piece, task$store, get("queue_model", envir = globalenv()))
}

## The queue: `workers` R processes (crew), each running one piece at a time;
## or, with workers = 0 (the dry run), every piece in this process, at once.
##   push(task)   adds a piece
##   busy()       whether pieces are still out
##   next_done()  waits for a piece to finish and returns its task with
##                its seconds; stops the run if the piece failed
##   stop()       ends the workers
start_queue <- function(workers) {
  # the pieces pushed and not yet handed back
  out <- new.env()
  out$tasks <- list()

  if (workers == 0) {
    return(list(
      push = function(task) {
        task$seconds <- run_task(task)
        out$tasks <- c(out$tasks, list(task))
      },
      busy = function() length(out$tasks) > 0,
      next_done = function() {
        task <- out$tasks[[1]]
        out$tasks <- out$tasks[-1]
        task
      },
      stop = function() invisible()
    ))
  }

  controller <- crew::crew_controller_local(
    workers = workers,
    # keep each worker's loaded code from one piece to the next
    reset_globals = FALSE,
    # a log per worker (its R output, including Stan's own messages), for
    # finding out why a piece of work failed
    options_local = crew::crew_options_local(log_directory = file.path("cache", "logs", "workers"))
  )
  controller$start()
  list(
    push = function(task) {
      name <- paste(task$design, piece_id(task$piece))
      out$tasks[[name]] <- task
      controller$push(
        command = run_task(task),
        data = list(task = task),
        globals = list(run_task = run_task),
        packages = c("flatcor", "dplyr"),
        name = name
      )
    },
    busy = function() length(out$tasks) > 0,
    next_done = function() {
      repeat {
        finished <- controller$pop()
        if (!is.null(finished)) break
        controller$wait(mode = "one", seconds_timeout = 60)
      }
      if (!identical(finished$status, "success")) {
        stop(
          "A piece of work failed (", finished$name, "): ", finished$error,
          "\nThe workers' logs are in cache/logs/workers/.",
          call. = FALSE
        )
      }
      task <- out$tasks[[finished$name]]
      out$tasks[[finished$name]] <- NULL
      task$seconds <- finished$result[[1]]
      task
    },
    stop = function() controller$terminate()
  )
}


## ---- Judging a cell ---------------------------------------------------------------

## The cells of a design: per cell, the rows of the plan it covers
## (the rows that share the values of `by`, as in next_plan())
cells_of <- function(plan, by) split(seq_len(nrow(plan)), do.call(paste, plan[by]))

## A table as the per-dataset file gives it back (write_per_dataset(), then
## read_per_dataset(), R/shared/per_dataset.R): a cell is judged on exactly
## the numbers its study's results are computed from. Text columns stay text:
## read alone, a cell's column of only "true" would come back as TRUE, where
## the study's whole file, with other values in it, keeps it as text (and a
## bootstrap's stream is named after those values).
as_read_back <- function(table) {
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file))
  table <- as.data.frame(table)
  table <- table[setdiff(names(table), c("variables", "chains"))]
  data.table::fwrite(table, file)
  text <- names(table)[vapply(table, is.character, logical(1))]
  as.data.frame(data.table::fread(
    file, na.strings = c("", "NA"), colClasses = list(character = text)
  ))
}

## The plan after judging one cell of a design: the cell's rows get the
## datasets it is estimated to lack (next_plan(), R/shared/rounds.R), or
## stay as they are when it has reached its targets. The cell is judged on
## its own results only: the tables of its study (design$table) computed
## from its own finished pieces.
judge_cell <- function(design, plan, rows, store, one_batch) {
  pieces <- pieces_of_plan(plan[rows, , drop = FALSE], design$batch_size)
  results <- bind_pieces(lapply(seq_len(nrow(pieces)), function(i) {
    read_piece(design$name, pieces[i, ], store)$rows
  }))
  columns <- c(design$by, "mcse", "target_mcse", "reached")
  tables <- design$table(as_read_back(results))
  table <- do.call(rbind, lapply(tables, function(t) as.data.frame(t)[columns]))
  next_plan(plan, table, design$by, design$batch_size, one_batch = one_batch)
}


## ---- The queue's run ----------------------------------------------------------------

## Runs every cell of `designs` (run_all.R's in_rounds) to its targets or to
## the budget. `workers`: as in start_queue(). `spent_elsewhere`: core hours
## the pipeline's own steps have cost so far (a run that continues).
## `one_extension`: at most one extra batch per cell (the smoke test, which
## tests the step, not the numbers).
run_rounds <- function(designs, store, budget_core_hours, workers, spent_elsewhere,
                       one_extension = FALSE) {
  queue <- start_queue(workers)
  on.exit(queue$stop())
  files <- code_files()

  # what changes while the queue runs, in one place
  run <- new.env()
  run$plans <- lapply(designs, function(design) {
    read_plan(design$name, store, design$cells, design$first_round)
  })
  # every piece's seconds, per design: those saved before (a run that
  # continues) and those that come back
  run$seconds <- lapply(designs, function(design) {
    saved <- list.files(file.path(store, "pieces", design$name), "[.]rds$", full.names = TRUE)
    vapply(saved, function(file) readRDS(file)$seconds, numeric(1), USE.NAMES = FALSE)
  })
  run$out <- list()        # per cell: the expected core hours of each piece out
  run$extended <- character(0)
  run$short <- list()      # per cell its next batch did not fit: what it would cost

  # core hours: every piece runs Stan's 4 chains, so 4 cores
  spent <- function() spent_elsewhere + sum(unlist(run$seconds)) * 4 / 3600
  out_cost <- function() sum(unlist(run$out))
  per_dataset <- function(name) {
    if (length(run$seconds[[name]]) == 0) return(0)
    mean(run$seconds[[name]]) * 4 / 3600 / designs[[name]]$batch_size
  }
  # the pieces of a cell's plan not run yet
  not_run <- function(name, rows) {
    pieces <- pieces_of_plan(run$plans[[name]][rows, , drop = FALSE], designs[[name]]$batch_size)
    missing <- vapply(seq_len(nrow(pieces)), function(i) {
      is.null(read_piece(name, pieces[i, ], store))
    }, logical(1))
    pieces[missing, , drop = FALSE]
  }
  push_pieces <- function(name, key, pieces) {
    cell <- paste(name, key)
    for (i in seq_len(nrow(pieces))) {
      run$out[[cell]] <- c(run$out[[cell]], (pieces$to[i] - pieces$from[i] + 1) * per_dataset(name))
      queue$push(list(design = name, key = key, piece = pieces[i, ], store = store, files = files))
    }
  }

  # a cell with no pieces out: run what its plan still lacks (a run that
  # continues), or judge it and start its next batch
  next_step <- function(name, key) {
    design <- designs[[name]]
    cell <- paste(name, key)
    rows <- cells_of(run$plans[[name]], design$by)[[key]]
    missing <- not_run(name, rows)
    if (nrow(missing) > 0) return(push_pieces(name, key, missing))

    plan_next <- judge_cell(design, run$plans[[name]], rows, store, one_batch = one_extension)
    added <- sum(plan_next$datasets[rows] - run$plans[[name]]$datasets[rows])
    if (added == 0 || (one_extension && cell %in% run$extended)) {
      datasets <- sum(run$plans[[name]]$datasets[rows])
      return(message(sprintf("%s: done, %d datasets", cell, datasets)))
    }
    cost <- added * per_dataset(name)
    if (spent() + out_cost() + cost > budget_core_hours) {
      run$short[[cell]] <- cost
      return(message(sprintf("%s: next batch (%.1f core hours) not started: budget", cell, cost)))
    }
    message(sprintf("%s: %d more datasets", cell, added))
    run$extended <- c(run$extended, cell)
    run$plans[[name]]$datasets[rows] <- plan_next$datasets[rows]
    write_plan(run$plans[[name]], name, store)
    push_pieces(name, key, not_run(name, rows))
  }

  # every cell starts; then each finished piece may complete its cell's batch
  for (name in names(designs)) {
    for (key in names(cells_of(run$plans[[name]], designs[[name]]$by))) next_step(name, key)
  }
  while (queue$busy()) {
    task <- queue$next_done()
    cell <- paste(task$design, task$key)
    run$seconds[[task$design]] <- c(run$seconds[[task$design]], task$seconds)
    run$out[[cell]] <- run$out[[cell]][-1]
    if (length(run$out[[cell]]) == 0) next_step(task$design, task$key)
  }

  if (length(run$short) == 0) {
    message(sprintf("Every cell is done. %.1f core hours so far.", spent()))
  } else {
    message(sprintf(
      paste(
        "%.1f core hours so far. %d cells did not start their next batch: together",
        "about %.1f core hours, more than the %.1f left of the budget of %.1f. To run",
        "them, raise the budget (VALIDATION_BUDGET_CORE_HOURS) and start again; to stop",
        "here, the cells short of a target are marked in their study's results (column",
        "reached)."
      ),
      spent(), length(run$short), sum(unlist(run$short)),
      max(0, budget_core_hours - spent()), budget_core_hours
    ))
  }
  invisible(run$plans)
}
