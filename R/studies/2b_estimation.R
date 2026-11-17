## =====================================================================
##  Study 2b: posterior estimation, does each prior behave like a default?
##  (docs/working/REDESIGN.md: heading 2, row b)
##
##  Question. On data with a known truth, how close are each prior's
##  posterior estimates to it: the bias of a slope's posterior mean, the
##  coverage of its 95% interval, and how well the posterior means predict
##  new data compared with ordinary least squares.
##  Data. Data set A (R/data/data_set_A.R), 16 cells: k = 2 and 10, n = 50
##  and 200, predictor correlation 0 and .7, the other slopes sparse or
##  dense; the slope of x1 is 0.2.
##  Methods. The three priors (flatcor's Stan model, fit_regression()), the
##  two references JZS (scale sqrt(2)/4) and Zellner-Siow (scale 1), both from
##  BayesFactor's posterior sampler, and OLS (lm()).
##  Estimands, per method and cell:
##    bias        (mean of x1's posterior mean - 0.2) / 0.2; reference 0 (OLS
##                is unbiased; a prior that shrinks pulls it below 0)
##    coverage    share of 95% intervals for x1 that contain 0.2; reference .95
##                (nominal; not guaranteed for a Bayesian interval when the
##                data come from fixed slopes rather than from the prior)
##    prediction  mean excess prediction error, (b - beta)' Sigma (b - beta),
##                divided by OLS's; reference 1 (below 1: better than OLS)
##  Reporting (docs/working/REDESIGN.md, the reporting rule). No pass or fail:
##  every estimate is shown with its reference, its Monte Carlo SE and its
##  distance from the reference in Monte Carlo SEs, and per quantity how many
##  lie beyond 2 and 3 MCSEs against what chance gives.
##  Stan's diagnostics (divergent transitions) are reported as counts.
##  n_sim. 475 datasets per cell (R/shared/precision.R, coverage_per_cell):
##  a coverage near .95 has MCSE <= .01. MCSEs: bias sd/sqrt(n_sim);
##  coverage sqrt(p (1 - p) / n_sim); the prediction ratio by a paired
##  bootstrap (Morris et al., 2019).
##
##  Run from the project folder:  Rscript R/studies/2b_estimation.R
##  Smoke test (one cell, 10 datasets):
##    VALIDATION_SMOKE=1 Rscript R/studies/2b_estimation.R
##  Raw results are cached in cache/study_2b/flatcor_<build>/ (delete it to
##  rerun); the tables are saved to results/studies/study_2b.rds. The Stan
##  fits run 4 chains in parallel, so SimDesign uses one worker per 4 cores.
## =====================================================================

library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, save_figure(), save_study_results()
source("R/shared/comparison_methods.R")   # contenders, method_scales()
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion()
source("R/data/data_sets.R")              # with_stream_of()
source("R/data/data_set_A.R")             # generate_A(), slopes_A(), cell_id_A()

study_folder <- study_cache_folder("2b")
bootstrap_replicates <- if (smoke_mode) 100 else 1000

design <- expand.grid(k = c(2, 10), n = c(50, 200), rho = c(0, 0.7),
                      others = c("sparse", "dense"), slope_x1 = 0.2, stringsAsFactors = FALSE)
if (smoke_mode) design <- subset(design, k == 2 & n == 50 & rho == 0 & others == "sparse")
n_sim <- n_sim_of("coverage_per_cell")
true_slope_x1 <- 0.2

## Stan: 4 chains of 1,000 draws after 1,000 warm-up; short chains in smoke mode
stan_chains   <- if (smoke_mode) 2 else 4
stan_warmup   <- if (smoke_mode) 200 else 1000
stan_sampling <- if (smoke_mode) 200 else 1000
bayesfactor_iterations <- if (smoke_mode) 1000 else 4000
workers <- max(1, number_of_cores %/% stan_chains)

methods <- c(names(contenders), "JZS default", "Zellner-Siow", "OLS")
method_key <- setNames(make.names(methods), methods)   # without spaces, for SimDesign's columns


## ---- The simulation ------------------------------------------------------------
analyse_dataset <- function(condition, dat, fixed_objects) {
  dataset_id  <- paste(cell_id_A(condition), condition$REPLICATION, sep = "-r")
  k <- condition$k
  covariance  <- matrix(condition$rho, k, k); diag(covariance) <- 1
  true_slopes <- slopes_A(k, condition$others, condition$slope_x1)
  columns <- colnames(dat$predictors)
  data    <- data.frame(dat$predictors, y = dat$outcome)
  formula <- reformulate(columns, response = "y")
  ## every fit's seed follows from its dataset's ID, so results do not depend
  ## on the number of cores or the order of the fits
  seed_of <- function(label) abs(digest::digest2int(paste(dataset_id, label))) %% .Machine$integer.max

  ## x1's estimate and 95% interval, and the excess prediction error of the
  ## estimated slopes
  summary_of <- function(estimates, lower, upper) {
    error <- estimates - true_slopes
    c(estimate = estimates[[1]], lower = lower, upper = upper,
      prediction_error = drop(t(error) %*% covariance %*% error))
  }
  from_draws <- function(draws)
    summary_of(colMeans(draws), quantile(draws[, 1], 0.025, names = FALSE),
               quantile(draws[, 1], 0.975, names = FALSE))

  per_method <- list(); divergent <- numeric(0)
  for (method in names(contenders)) {
    fit <- suppressWarnings(fit_regression(formula, data, prior = contenders[[method]],
                                           chains = stan_chains, iter_warmup = stan_warmup,
                                           iter_sampling = stan_sampling, seed = seed_of(method)))
    divergent[method_key[[method]]] <- fit$divergent
    per_method[[method]] <- from_draws(as.matrix(as.data.frame(fit$draws)[, columns, drop = FALSE]))
  }
  bayesfactor_draws <- function(scale, label)
    with_stream_of(paste(dataset_id, label), as.matrix(as.data.frame(
      BayesFactor::lmBF(formula, data, rscaleCont = scale, posterior = TRUE,
                        iterations = bayesfactor_iterations, progress = FALSE))[, columns, drop = FALSE]))
  per_method[["JZS default"]]  <- from_draws(bayesfactor_draws(sqrt(2) / 4, "JZS"))
  per_method[["Zellner-Siow"]] <- from_draws(bayesfactor_draws(1, "ZS"))
  ols <- lm(formula, data)
  per_method[["OLS"]] <- summary_of(coef(ols)[columns], confint(ols)["x1", 1], confint(ols)["x1", 2])

  c(unlist(lapply(methods, function(method)
      setNames(per_method[[method]], paste(names(per_method[[method]]), method_key[[method]], sep = ".")))),
    setNames(divergent, paste0("divergent.", names(divergent))))
}

## SimDesign needs a summarise function to keep the per-dataset results
summarise_cell <- function(condition, results, fixed_objects) c(datasets = nrow(results))

simulation_file <- file.path(study_folder, "estimation.rds")
simulation <- if (file.exists(simulation_file)) {
  message("Using cached results: ", simulation_file); readRDS(simulation_file)
} else {
  result <- runSimulation(design = design, replications = n_sim,
    generate = function(condition, fixed_objects) generate_A(condition, condition$REPLICATION),
    analyse = analyse_dataset, summarise = summarise_cell,
    parallel = workers > 1, ncores = workers, packages = c("flatcor", "BayesFactor"),
    control = list(include_replication_index = TRUE), store_results = TRUE, progress = FALSE,
    save_details = list(out_rootdir = study_folder, tmpfilename = "resume_estimation.qs2"))
  saveRDS(result, simulation_file)
  result
}

## one row per dataset and method
per_dataset <- SimExtract(simulation, what = "results")
estimates <- per_dataset |>
  select(-starts_with("divergent.")) |>
  pivot_longer(matches("^(estimate|lower|upper|prediction_error)[.]"),
               names_to = c(".value", "key"), names_pattern = "^([a-z_]+)[.](.*)$") |>
  mutate(method = names(method_key)[match(key, method_key)],
         method = factor(method, levels = methods),
         cell = sprintf("k = %d, n = %d, rho = %.1f, %s", k, n, rho, others))


## ---- Bias and coverage of x1 ---------------------------------------------------
bias_and_coverage <- estimates |>
  group_by(k, n, rho, others, cell, method) |>
  dplyr::summarise(datasets = n(),
                   bias = (mean(estimate) - true_slope_x1) / true_slope_x1,
                   bias_mcse = sd(estimate) / sqrt(datasets) / true_slope_x1,
                   coverage = mean(lower <= true_slope_x1 & true_slope_x1 <= upper),
                   coverage_mcse = mcse_of_proportion(coverage, datasets), .groups = "drop")


## ---- Excess prediction error relative to OLS, with a paired bootstrap SE --------
RNGkind("L'Ecuyer-CMRG"); set.seed(9)   # the bootstrap's own, fixed stream
prediction <- estimates |>
  group_by(k, n, rho, others, cell) |>
  group_modify(function(rows, keys) {
    error_of <- function(name)                     # in replication order
      rows$prediction_error[rows$method == name][order(rows$REPLICATION[rows$method == name])]
    ols <- error_of("OLS")
    resamples <- replicate(bootstrap_replicates, sample(length(ols), replace = TRUE), simplify = FALSE)
    bind_rows(lapply(setdiff(methods, "OLS"), function(name) {
      errors <- error_of(name)
      tibble(method = name, ratio = mean(errors) / mean(ols),
             ratio_mcse = sd(vapply(resamples, function(r) mean(errors[r]) / mean(ols[r]), numeric(1))))
    }))
  }) |>
  ungroup()


## ---- One table: estimate, reference, Monte Carlo SE, distance in MCSEs ---------
results_table <- bind_rows(
  transmute(bias_and_coverage, k, n, rho, others, cell, method, quantity = "bias of x1 (share of 0.2)",
            reference = 0, estimate = bias, mcse = bias_mcse),
  transmute(bias_and_coverage, k, n, rho, others, cell, method, quantity = "coverage of x1's 95% interval",
            reference = 0.95, estimate = coverage, mcse = coverage_mcse),
  transmute(prediction, k, n, rho, others, cell, method, quantity = "prediction error / OLS",
            reference = 1, estimate = ratio, mcse = ratio_mcse)) |>
  mutate(distance_in_mcse = (estimate - reference) / mcse)

divergent <- per_dataset |>
  select(starts_with("divergent.")) |>
  pivot_longer(everything(), names_to = "key", names_prefix = "divergent.", values_to = "divergent") |>
  mutate(method = names(method_key)[match(key, method_key)]) |>
  group_by(method) |>
  dplyr::summarise(fits = n(), fits_with_divergences = sum(divergent > 0),
                   divergent_transitions = sum(divergent), .groups = "drop")

for (q in unique(results_table$quantity)) {
  cat(sprintf("\n%s (reference %g): estimate (MCSE)\n", q, results_table$reference[results_table$quantity == q][1]))
  print(as.data.frame(results_table |> filter(quantity == q) |>
    mutate(value = sprintf("%.3f (%.3f)", estimate, mcse)) |>
    select(cell, method, value) |>
    pivot_wider(names_from = method, values_from = value)), row.names = FALSE)
}
## how many estimates lie beyond 2 and 3 Monte Carlo SEs of their reference,
## against what chance alone gives if the reference held (about 5% and 0.3%);
## here the references are not guaranteed, so this describes, not tests
beyond_reference <- results_table |>
  group_by(quantity, method) |>
  dplyr::summarise(estimates = n(), beyond_2_mcse = sum(abs(distance_in_mcse) > 2),
                   beyond_3_mcse = sum(abs(distance_in_mcse) > 3), .groups = "drop") |>
  mutate(expected_by_chance_2 = round(0.0455 * estimates, 1),
         expected_by_chance_3 = round(0.0027 * estimates, 1))
cat("\nEstimates beyond 2 and 3 Monte Carlo SEs of their reference:\n")
print(as.data.frame(beyond_reference), row.names = FALSE)
cat("\nStan's divergent transitions:\n"); print(as.data.frame(divergent), row.names = FALSE)


## ---- Figures -------------------------------------------------------------------
cells_facet <- facet_grid(others + rho ~ k,
  labeller = labeller(k = function(v) paste("k =", v), rho = function(v) paste("correlation =", v)))
figure_of <- function(q, ylab, title, methods_shown = methods)
  ggplot(filter(results_table, quantity == q, method %in% methods_shown), aes(n, estimate, colour = method)) +
    geom_hline(aes(yintercept = reference), linetype = "dotted") +
    geom_pointrange(aes(ymin = estimate - 2 * mcse, ymax = estimate + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) + method_scales(methods_shown) + cells_facet +
    labs(x = "sample size (log scale)", y = ylab, title = title,
         subtitle = sprintf("%d datasets per cell; bars: 2 Monte Carlo SE; dotted: the reference", n_sim))
save_figure(figure_of("coverage of x1's 95% interval", "coverage of x1's 95% interval",
                      "Do the intervals contain the true slope of x1 (0.2)?"),
            "study_2b_coverage", width = 9, height = 10)
save_figure(figure_of("bias of x1 (share of 0.2)", "bias of x1, as a share of its true value",
                      "Bias: how far estimates of x1 are pulled towards zero"),
            "study_2b_bias", width = 9, height = 10)
save_figure(figure_of("prediction error / OLS", "excess prediction error / OLS",
                      "Do the posterior means predict new data better than OLS?", setdiff(methods, "OLS")),
            "study_2b_prediction", width = 9, height = 10)

save_study_results("2b", list(estimates = results_table, beyond_reference = beyond_reference,
                              divergent = divergent))
