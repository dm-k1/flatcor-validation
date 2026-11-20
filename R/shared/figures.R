## =====================================================================
##  Shared by the figure scripts in R/3_figures/: where a study's results file
##  and per-dataset files are, and where its figures go (cache/ for the smoke
##  test and the dry run)
## =====================================================================

smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
dry <- identical(Sys.getenv("VALIDATION_DRY"), "1")

## The folder of the results: results/ (cache/ for the smoke test and the dry run)
results_root <- if (dry) {
  file.path("cache", "dry_results")
} else if (smoke) {
  file.path("cache", "smoke_results")
} else {
  "results"
}

## The results file of a study (results/studies/study_<cell>.rds)
read_results <- function(study) {
  readRDS(file.path(results_root, "studies", sprintf("study_%s.rds", study)))
}

## A per-dataset file of a study (results/per_dataset/study_<file>.csv or
## .csv.gz, R/shared/per_dataset.R)
read_per_dataset_file <- function(file) {
  path <- list.files(
    file.path(results_root, "per_dataset"),
    pattern = sprintf("^study_%s[.]csv([.]gz)?$", file), full.names = TRUE
  )
  stopifnot(length(path) == 1)
  as.data.frame(data.table::fread(path, na.strings = c("", "NA")))
}

## Where the figures go
figure_folder <- if (dry) {
  file.path("cache", "dry_figures")
} else if (smoke) {
  file.path("cache", "smoke_figures")
} else {
  file.path("results", "figures")
}
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)
