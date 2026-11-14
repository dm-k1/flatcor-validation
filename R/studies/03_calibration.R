## =====================================================================
##  Study 3: are the posterior and the Bayes factor calibrated?
##  (docs/appendix/C_study_design.md, study 3)
##
##  Aims. Two calibration checks on data drawn from each prior, the same for
##  every prior:
##    Part A, the posterior: simulation-based calibration (Talts et al., 2018)
##    of flatcor's Stan model, with the SBC package. If the sampler is right,
##    the rank of each true value among its posterior draws is uniform.
##    Part B, the Bayes factor (Schad et al., 2023): among datasets with a
##    given BF10, the share drawn from the model without x1 equals the
##    posterior probability of that model, 1 / (1 + BF10), with equal prior
##    odds.
##  Data. Set B (R/data/data_set_B.R): predictors with correlation .3, the
##  intercept from N(0, 10) and the noise SD from lognormal(0, 1) (the Stan
##  model's proper stand-ins), the coefficients from the generating prior.
##    Part A: k = 8, n = 60 and k = 3, n = 50; all coefficients drawn.
##    Part B: k = 3 and 8, n = 50 and 200; each of the three priors
##    generates. Odd replications come from the model without x1 (x1 = 0, the
##    others drawn from that model's own prior for k - 1 predictors), even
##    ones from the model with all k (every coefficient drawn).
##  Estimands. Part A: the ranks of the true R-squared, b[1] and sigma, and
##  the divergent transitions of each fit. Part B:
##  the observed share of datasets from the model without x1, per bin of log
##  BF10, against the mean predicted share; the logistic calibration line of
##  "without x1" on log BF01 (intercept 0, slope 1 when calibrated).
##  Methods. Part A: flatcor's Stan model under each prior. Part B: the Bayes
##  factor of the model with x1 against the model without it,
##  bayes_factor_models()'s route (log_bayes_factor_vs_intercept_only() of
##  each model), under the prior that generated the data.
##  Performance. Part A: the SBC package's ranks (thinned so each fit gives
##  100) and ECDF-difference plots; uniformity per quantity by a chi-square
##  over 20 rank bins. Part B: per bin, the observed share with its Monte
##  Carlo SD under calibration; a chi-square over the bins; the logistic slope
##  with its SE.
##  n_sim. Part A: 1,000 fits per prior and setting (precision.R, sbc_fits).
##  Part B: 20,000 datasets per cell (precision.R, bf_calibration): bins
##  holding at least 625 datasets have shares with MCSE <= .02.
##
##  Checks:
##    3.1 Under each prior and setting, the ranks of R-squared, b[1] and sigma
##        are uniform: chi-square over 20 bins, p > .01.
##    3.2 Under each prior and setting, the Stan model samples without
##        divergent transitions (SBC's backend diagnostics).
##    3.3 Each prior's Bayes factor is calibrated on data from its own prior:
##        the chi-square over the bins does not reject (p > .01) and the
##        logistic slope is within 3 SE of 1, in every cell.
##
##  Run from the project folder:  Rscript R/studies/03_calibration.R
##  Smoke test: VALIDATION_SMOKE=1 Rscript R/studies/03_calibration.R
##  Raw results are cached in cache/study_03/flatcor_<build>/ (delete it to
##  rerun); the tables and checks are saved to results/studies/study_03.rds.
## =====================================================================
library(SBC)
library(SimDesign)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # contenders
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion()
source("R/data/data_sets.R")              # with_stream_of()
source("R/data/data_set_B.R")             # generate_B()

study_folder <- study_cache_folder("03")
future::plan(future::multisession, workers = number_of_cores)


## =====================================================================
## PART A: THE POSTERIOR (simulation-based calibration)
## =====================================================================
n_fits <- n_sim_of("sbc_fits")
settings <- data.frame(k = c(8, 3), n = c(60, 50))
if (smoke_mode) settings <- settings[1, , drop = FALSE]
iterations <- if (smoke_mode) c(warmup = 300, sampling = 300) else c(warmup = 1000, sampling = 1000)

## the executable goes to this study's cache, never into the installed package
## (read-only in an Apptainer image)
direction_model <- cmdstanr::cmdstan_model(system.file("stan", "flatcor.stan", package = "flatcor"),
                                           dir = study_folder)

## One SBC dataset: the true values (R-squared, b[1], sigma) and the Stan data
sbc_dataset <- function(prior, k, n, replication) {
  dat <- generate_B(prior, k, n, x1 = "prior", purpose = "sbc", replication = replication)
  list(variables = list(R_squared = dat$r_squared, `b[1]` = dat$coefficients[1], sigma = dat$noise_sd),
       generated = list(N = n, K = k, X = dat$predictors, y = dat$outcome, prior_only = 0L,
                        proper_stand_ins = 1L, r_squared_shape = prior_shape(prior, k)))
}

run_sbc <- function(model, prior, k, n, label) {
  cache_file <- file.path(study_folder, sprintf("sbc_run_%s_k%d_n%d.rds", label, k, n))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  datasets <- do.call(bind_datasets, lapply(seq_len(n_fits), function(replication) {
    one <- sbc_dataset(prior, k, n, replication)
    SBC_datasets(variables = posterior::as_draws_matrix(posterior::as_draws_rvars(
                   lapply(one$variables, posterior::rvar))), generated = list(one$generated))
  }))
  backend <- SBC_backend_cmdstan_sample(model, chains = 1, iter_warmup = iterations[["warmup"]],
                                        iter_sampling = iterations[["sampling"]], seed = 3,
                                        refresh = 0, show_messages = FALSE)
  result <- compute_SBC(datasets, backend, keep_fits = FALSE,
                        thin_ranks = iterations[["sampling"]] / 100)
  run <- list(stats = result$stats |> mutate(prior = label, k = k, n = n),
              diagnostics = as.data.frame(result$backend_diagnostics) |> mutate(prior = label, k = k, n = n))
  saveRDS(run, cache_file)
  run
}

runs <- unlist(lapply(seq_len(nrow(settings)), function(row) lapply(unname(contenders),
  function(prior) run_sbc(direction_model, prior, settings$k[row], settings$n[row], prior))),
  recursive = FALSE)
stats <- bind_rows(lapply(runs, `[[`, "stats"))
divergences <- bind_rows(lapply(runs, `[[`, "diagnostics")) |>
  group_by(prior, k, n) |>
  summarise(fits = n(), fits_with_divergences = sum(n_divergent > 0),
            divergent_transitions = sum(n_divergent), .groups = "drop")


## ---- Uniformity of the ranks, and the checks ----------------------------------------
uniformity <- stats |>
  group_by(prior, k, n, variable) |>
  summarise(fits = n(), p_value = {
              bins <- table(cut(rank / max_rank, breaks = seq(0, 1, length.out = 21), include.lowest = TRUE))
              suppressWarnings(chisq.test(bins)$p.value)
            }, .groups = "drop")
cat("\nRank uniformity (chi-square over 20 bins):\n")
print(as.data.frame(uniformity), digits = 3, row.names = FALSE)

for (row in seq_len(nrow(uniformity))) with(uniformity[row, ],
  record_check(sprintf("3.1 %s, k = %d, n = %d: ranks of %s uniform", prior, k, n, variable),
               p_value > 0.01, sprintf("chi-square p = %.3f over %d fits", p_value, fits)))

cat("\nDivergent transitions of the SBC fits:\n")
print(as.data.frame(divergences), row.names = FALSE)
for (row in seq_len(nrow(divergences))) with(divergences[row, ],
  record_check(sprintf("3.2 %s, k = %d, n = %d: no divergent transitions", prior, k, n),
               divergent_transitions == 0,
               sprintf("%d divergent transitions in %d of %d fits",
                       divergent_transitions, fits_with_divergences, fits)))

## one panel row per prior (plot_ecdf_diff() needs one set of datasets at a time)
save_figure(
  patchwork::wrap_plots(lapply(unique(stats$prior), function(prior_name)
    plot_ecdf_diff(as.data.frame(filter(stats, prior == prior_name, k == settings$k[1]))) +
      labs(title = prior_name)), ncol = 1) +
    patchwork::plot_annotation(title = sprintf(
      "Simulation-based calibration, k = %d, n = %d: ECDF difference of the ranks",
      settings$k[1], settings$n[1])),
  "study_03_sbc_ecdf_difference", width = 10, height = 14)


## =====================================================================
## PART B: THE BAYES FACTOR
## =====================================================================
design <- createDesign(generating_prior = unname(contenders), k = c(3, 8), n = c(50, 200))
if (smoke_mode) design <- subset(design, k == 3 & n == 50)
n_sim <- if (smoke_mode) 200L else n_sim_of("bf_calibration")   # calibration needs a few hundred to bin


## ---- The simulation ------------------------------------------------------------
## The Bayes factor of the model with x1 against the model without it: each
## model against the intercept-only model, divided (bayes_factor_models()'s
## route)
simulate_data <- function(condition, fixed_objects) {
  r <- condition$REPLICATION
  generate_B(condition$generating_prior, condition$k, condition$n,
             x1 = if (r %% 2 == 1) "reduced" else "prior", purpose = "calibration", replication = r)
}

## log BF10 of the model with x1 against the model without it, under the prior
## that generated the data
analyse_dataset <- function(condition, dat, fixed_objects) {
  n <- condition$n; k <- condition$k; prior <- condition$generating_prior
  r_squared_with    <- sample_r_squared(dat$predictors, dat$outcome)
  r_squared_without <- sample_r_squared(dat$predictors[, -1, drop = FALSE], dat$outcome)
  c(is_null = as.numeric(dat$coefficients[1] == 0),
    log_BF10 = log_bayes_factor_vs_intercept_only(r_squared_with, n, k, prior) -
      log_bayes_factor_vs_intercept_only(r_squared_without, n, k - 1, prior))
}

summarise_cell <- function(condition, results, fixed_objects) c(share_null = mean(results$is_null))

results_file <- file.path(study_folder, "bayes_factor_calibration.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset, summarise = summarise_cell,
                              parallel = number_of_cores > 1, ncores = number_of_cores,
                              packages = c("flatcor"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE,
                              save_details = list(out_rootdir = study_folder, tmpfilename = "resume_bayes_factor_calibration.qs2"))
  saveRDS(simulation, results_file)
}

## one row per dataset, tested under the prior that generated it
calibration_runs <- SimExtract(simulation, what = "results") |>
  rename(prior = generating_prior) |>
  mutate(is_null = is_null == 1,
         predicted_null = plogis(-log_BF10),          # 1 / (1 + BF10)
         cell = sprintf("k = %d, n = %d", k, n))


## ---- Per bin: observed share of nulls against the mean of 1 / (1 + BF10) --------
bin_edges <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, Inf)
calibration_bins <- calibration_runs |>
  mutate(bin = cut(log_BF10, bin_edges)) |>
  group_by(prior, cell, bin) |>
  summarise(datasets = n(), observed = mean(is_null), predicted = mean(predicted_null),
            ## Monte Carlo SD of the observed share if calibrated (Bernoulli(p_i) draws)
            mc_sd = sqrt(sum(predicted_null * (1 - predicted_null))) / n(),
            mcse_observed = mcse_of_proportion(mean(is_null), n()),
            .groups = "drop") |>
  mutate(z = (observed - predicted) / mc_sd)

## Per cell: the chi-square over the bins with at least 5 expected nulls and 5
## expected effects, and the logistic calibration line
chi_square_tests <- calibration_bins |>
  filter(datasets * predicted >= 5, datasets * (1 - predicted) >= 5) |>
  group_by(prior, cell) |>
  summarise(bins = n(), chi_square = sum(z^2), largest_abs_z = max(abs(z)), .groups = "drop") |>
  mutate(p_value = pchisq(chi_square, bins, lower.tail = FALSE))
logistic_lines <- calibration_runs |>
  group_by(prior, cell) |>
  group_modify(function(rows, key) {
    ## (datasets with huge Bayes factors give fitted probabilities of 0 or 1: harmless)
    fit <- suppressWarnings(glm(is_null ~ I(-log_BF10), family = binomial, data = rows))
    estimates <- summary(fit)$coefficients
    tibble(datasets = nrow(rows), intercept = estimates[1, 1], intercept_se = estimates[1, 2],
           slope = estimates[2, 1], slope_se = estimates[2, 2])
  }) |>
  ungroup()
calibration_tests <- left_join(chi_square_tests, logistic_lines, by = c("prior", "cell"))

cat("\nCalibration per bin (share of nulls: observed with its MCSE, predicted = mean 1/(1 + BF10)):\n")
print(as.data.frame(calibration_bins |> select(prior, cell, bin, datasets,
                                               observed, mcse_observed, predicted)), digits = 3)
cat("\nChi-square over the bins and the logistic calibration line per cell:\n")
print(as.data.frame(calibration_tests), digits = 3)


## ---- Checks --------------------------------------------------------------------
for (prior_name in unname(contenders)) for (cell_name in unique(calibration_tests$cell)) {
  row <- filter(calibration_tests, prior == prior_name, cell == cell_name)
  record_check(sprintf("3.3 %s: Bayes factor calibrated on its own prior, %s", prior_name, cell_name),
               row$p_value > 0.01 && abs(row$slope - 1) < 3 * row$slope_se,
               sprintf("chi2 = %.1f on %d bins (p = %.3f; %d datasets); slope %.3f (SE %.3f)",
                       row$chi_square, row$bins, row$p_value, row$datasets, row$slope, row$slope_se))
}


## ---- Figure --------------------------------------------------------------------
save_figure(
  ggplot(filter(calibration_bins, datasets >= 20), aes(predicted, observed)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = predicted - 2 * mc_sd, ymax = predicted + 2 * mc_sd),
                  width = 0, colour = "grey70") +
    geom_point(size = 2, colour = "#0072B2") + geom_line(colour = "#0072B2") +
    facet_grid(prior ~ cell) +
    coord_equal() +
    labs(x = "predicted share without x1, mean of 1 / (1 + BF10)", y = "observed share without x1",
         title = "Calibration of the Bayes factor for x1 (bins of log BF10, grey: 2 Monte Carlo SD)",
         subtitle = sprintf("%d datasets per cell; rows: the prior that generated the data and computes the Bayes factor", n_sim)),
  "study_03_bayes_factor_calibration", width = 12, height = 10)

save_study_results("03", list(sbc_stats = stats, sbc_uniformity = uniformity, sbc_divergences = divergences,
                              calibration_bins = calibration_bins, calibration_tests = calibration_tests))
report_checks()
