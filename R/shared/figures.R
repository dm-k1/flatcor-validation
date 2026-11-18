## =====================================================================
##  Shared by the figure scripts in R/3_figures/: where a study's results file
##  is, and where its figures go (cache/ for the smoke test)
## =====================================================================

smoke <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")

## The results file of a study (results/studies/study_<cell>.rds)
read_results <- function(study) {
  folder <- if (smoke) file.path("cache", "smoke_results") else "results"
  readRDS(file.path(folder, "studies", sprintf("study_%s.rds", study)))
}

## Where the figures go
figure_folder <- if (smoke) file.path("cache", "smoke_figures") else file.path("results", "figures")
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)
