## =====================================================================
##  The per-dataset results of the final validation (v1), as CSV files
##
##  Run once on 20 Nov 2026 from copies of the two final targets stores
##  (S3: real-1a2b/store for 1a, 2a, 2b; real-3a3b/store for 1b, 3a, 3b),
##  from the repository root:
##    Rscript results/per_dataset/export_v1.R <real-1a2b store> <real-3a3b store>
##  It reads what the runs stored and writes it out unchanged: no step of the
##  pipeline runs and no number is computed. Where a results file in
##  results/studies/ already holds the per-dataset table, that table is
##  written; 1a's prior draws and 3b's log Bayes factors come from the
##  stores. The checks at the end compare the stored rows with the
##  results files they were summarised into. See README.md here.
## =====================================================================

stores <- commandArgs(trailingOnly = TRUE)
if (length(stores) != 2) stop("Give the two stores: <real-1a2b store> <real-3a3b store>")
store_1a2b <- stores[1]
store_3a3b <- stores[2]
folder <- file.path("results", "per_dataset")

from_store <- function(name, store) targets::tar_read_raw(name, store = store)
results_table <- function(study, table) {
  readRDS(file.path("results", "studies", sprintf("study_%s.rds", study)))$tables[[table]]
}
write_table <- function(table, file) {
  path <- file.path(folder, file)
  utils::write.csv(as.data.frame(table), path, row.names = FALSE)
  cat(sprintf("%-40s %7d rows\n", file, nrow(table)))
}

## ---- 1a: the prior draws, one row per draw ----------------------------------------------
samples <- from_store("samples_1a", store_1a2b)
draws_1a <- do.call(rbind, lapply(samples, function(setting) {
  do.call(rbind, lapply(names(setting$quantities), function(quantity) {
    draws <- setting$quantities[[quantity]]
    data.frame(
      prior = setting$setting$prior,
      design = setting$setting$design,
      quantity = quantity,
      chain = rep(seq_len(ncol(draws)), each = nrow(draws)),
      iteration = rep(seq_len(nrow(draws)), times = ncol(draws)),
      value = as.vector(draws)
    )
  }))
}))
draws_file <- gzfile(file.path(folder, "study_1a_prior_draws.csv.gz"))
utils::write.csv(draws_1a, draws_file, row.names = FALSE)
cat(sprintf("%-40s %7d rows\n", "study_1a_prior_draws.csv.gz", nrow(draws_1a)))

## ---- 1b, 2a, 2b, 3a: the per-dataset tables of the results files ------------------------
write_table(results_table("1b", "invariance"), "study_1b_changes.csv")
write_table(results_table("2a", "sbc_stats"), "study_2a_sbc_fits.csv")
write_table(results_table("2b", "per_dataset"), "study_2b_fits.csv")
write_table(results_table("3a", "bayes_factors"), "study_3a_calibration.csv")
write_table(results_table("3a", "repeats"), "study_3a_repeated_runs.csv")

## ---- 3b: the log Bayes factors, from the store --------------------------------------------
formula_3b <- from_store("formula_3b", store_3a3b)
rates_3b <- from_store("rates_3b", store_3a3b)
growing_3b <- from_store("growing_3b", store_3a3b)
write_table(formula_3b, "study_3b_formula.csv")
write_table(rates_3b, "study_3b_rates.csv")
write_table(growing_3b, "study_3b_growing_k.csv")

## ---- Checks: the stored rows are the ones the results files summarise ---------------------
same <- function(a, b) {
  isTRUE(all.equal(as.data.frame(a), as.data.frame(b), check.attributes = FALSE))
}
stopifnot(
  "2b: the stored fits are the results file's per-dataset table" =
    same(from_store("fits_2b", store_1a2b), results_table("2b", "per_dataset")),
  "3a: the stored calibration runs are the results file's table" =
    same(from_store("calibration_runs_3a", store_3a3b), results_table("3a", "bayes_factors")),
  "3a: the stored repeated runs are the results file's table" =
    same(from_store("repeated_3a", store_3a3b), results_table("3a", "repeats"))
)
medians <- stats::aggregate(
  log_bf ~ truth + growth + n + method,
  data = growing_3b, FUN = stats::median
)
reported <- merge(results_table("3b", "growing_medians"), medians)
stopifnot(
  "3b: every reported median has its stored log BFs" =
    nrow(reported) == nrow(results_table("3b", "growing_medians")),
  "3b: the medians of the stored log BFs are the reported medians" =
    isTRUE(all.equal(reported$median_log_bf, reported$log_bf))
)
cat("Checks: the stored rows match the results files.\n")
