## =====================================================================
##  The per-dataset files: results/per_dataset/study_<cell>_<what>.csv(.gz)
##  (cache/smoke_results/per_dataset/ for the smoke test, cache/dry_results/
##  for the dry run, cache/timing/ for a timing run)
##
##  Three layers. Layer 1: every study writes, per dataset (and prior or
##  method), everything it computed, in one step per file (format = "file"):
##  the numbers, both log marginal likelihoods and bridge sampling's errors,
##  the summary of every Stan variable (file _variables), the sampler per
##  chain (file _chains), seeds, names and seconds. Layer 2: the metrics of
##  docs/validation_tests.md are computed from those files, read back here,
##  so the published files are exactly what the results come from. Layer 3:
##  the figures and the notebook read both.
##
##  A row's `id` is the name its random numbers come from (seed = a hash of
##  it, R/shared/random_streams.R): the dataset can be made again from it.
## =====================================================================

## The folder of the per-dataset files
per_dataset_folder <- function() {
  if (identical(Sys.getenv("VALIDATION_DRY"), "1")) {
    return(file.path("cache", "dry_results", "per_dataset"))
  }
  smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
  timing <- !is.na(suppressWarnings(as.numeric(Sys.getenv("VALIDATION_TIMING", NA))))
  root <- if (smoke) {
    file.path("cache", "smoke_results")
  } else if (timing) {
    file.path("cache", "timing", "results")
  } else {
    "results"
  }
  file.path(root, "per_dataset")
}

## The name a seed comes from (seed_of() hashes the same parts)
name_of <- function(...) paste(...)

## The machine a piece of work ran on, for its rows
machine_of <- function() paste(Sys.info()[["nodename"]], Sys.info()[["machine"]])

## Writes one study's table: its plain columns to study_<file>.csv (.csv.gz
## with gzip = TRUE, for draws), and, if
## it has them, its list columns `variables` and `chains` (a data frame per
## row, joined back by the row's `id`) to study_<file>_variables.csv.gz and
## study_<file>_chains.csv.gz. Returns the paths (a targets "file" step).
##   split_by  a column: those two files are written once per value of it
##             (study_<file>_variables_<value>.csv.gz), to keep each under
##             GitHub's 100 MB limit (3b's growing k: one per growth rule)
write_per_dataset <- function(table, file, gzip = FALSE, split_by = NULL) {
  folder <- per_dataset_folder()
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  table <- as.data.frame(table)
  # two columns of the same name (a design column and one of the Bayes
  # factor's) come out of bind_cols() renamed as name...3: a mistake
  clashes <- grep("[.][.][.][0-9]+$", names(table), value = TRUE)
  if (length(clashes) > 0) stop("Columns named twice in ", file, ": ", toString(clashes))
  nested <- intersect(c("variables", "chains"), names(table))
  paths <- file.path(folder, sprintf(if (gzip) "study_%s.csv.gz" else "study_%s.csv", file))
  data.table::fwrite(table[setdiff(names(table), nested)], paths)

  for (column in nested) {
    keep <- !vapply(table[[column]], is.null, logical(1))
    groups <- if (is.null(split_by)) "" else sort(unique(table[[split_by]]))
    for (group in groups) {
      rows <- if (group == "") keep else keep & table[[split_by]] == group
      long <- data.table::rbindlist(
        Map(
          function(id, part) cbind(id = id, as.data.frame(part)),
          table$id[rows], table[[column]][rows]
        ),
        fill = TRUE
      )
      # "k = sqrt(n)" -> "_k_sqrt_n"
      suffix <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", group))
      if (group != "") suffix <- paste0("_", suffix)
      path <- file.path(folder, sprintf("study_%s_%s%s.csv.gz", file, column, suffix))
      data.table::fwrite(long, path)
      paths <- c(paths, path)
    }
  }
  paths
}

## Reads a file written by write_per_dataset() back, as a data frame. `paths`:
## the step's paths; `part`: "" for the main file, or "variables", "chains"
## (the files of a split_by, bound back into one table)
read_per_dataset <- function(paths, part = "") {
  nested <- "_(variables|chains)(_[A-Za-z0-9_]+)?[.]csv[.]gz$"
  path <- if (part == "") {
    paths[!grepl(nested, paths)]
  } else {
    paths[grepl(sprintf("_%s(_[A-Za-z0-9_]+)?[.]csv[.]gz$", part), paths)]
  }
  stopifnot(length(path) >= 1, part != "" || length(path) == 1)
  tables <- lapply(path, data.table::fread, na.strings = c("", "NA"))
  as.data.frame(data.table::rbindlist(tables, fill = TRUE))
}

## The pieces of a study as one table (targets' branches, as a list or
## already bound), keeping the list columns
bind_pieces <- function(pieces) {
  if (is.data.frame(pieces)) return(pieces)
  dplyr::bind_rows(pieces)
}
