## =====================================================================
##  Study 5: how does each prior behave on data with a known truth?
##  (docs/appendix/C_study_design.md, study 5)
##
##  Aims. The operating characteristics of the three priors next to the two
##  references, the same for every method:
##    Part A, one coefficient: how often the Bayes factor for x1 (the model
##    with x1 against the model without it) supports the model without it
##    (BF01 > 3), supports x1 (BF10 > 3) or stays undecided; how that depends
##    on the number of predictors k; how well it separates a null from a
##    small effect.
##    Part B, a factor: the same for a factor with q + 1 levels (the model
##    with the factor against the model without it).
##    Part C, estimation: bias, interval coverage and prediction error of the
##    posterior.
##  Data. Part A: set A (R/data/data_sets.R), 72 designs (k, n, predictor
##  correlation, the other slopes) at slopes of x1 of 0, 0.1, 0.2 and 0.3.
##  Part B: set A extended with a factor (R/data/data_set_A_factors.R): q = 2
##  to 5, n = 50, 200, 1,000, 2 or 10 other predictors that explain
##  R-squared = .5, factor effect Cohen's f^2 = 0, .02, .15. Part C: 16 cells
##  of set A (k 2 and 10, n 50 and 200, correlation 0 and .7, sparse and
##  dense) at slope 0.2, the same datasets as part A's for those cells.
##  Estimands. Parts A and B, per method and cell: the share of each answer
##  at BF 3; the AUC for a null against a small effect, and its difference
##  from JZS. Part C: the relative bias of x1's posterior mean, the coverage
##  of its 95% interval, the excess prediction error relative to OLS.
##  Methods. The three priors (flatcor: Bayes factors from
##  log_bayes_factor_vs_intercept_only(), posteriors from fit_regression()),
##  JZS (BayesFactor, scale sqrt(2)/4; for a factor lmBF with its default
##  scales) and Zellner-Siow (scale 1); in part C also OLS. Every Bayes factor
##  compares two models, each with its own prior.
##  Performance. Shares with their Monte Carlo SE sqrt(p (1 - p) / n_sim);
##  AUC differences and prediction-error ratios with paired bootstrap SEs;
##  bias and coverage with rsimsum's MCSE.
##  n_sim. Parts A and B: 625 datasets per cell (precision.R, rate_per_cell):
##  every per-cell share has MCSE <= .02. Part C: 475 per cell
##  (coverage_per_cell): coverage near .95 has MCSE <= .01.
##
##  Checks, for every method (the same question for each):
##    5.1 Part A: the share of true nulls the method confirms is the same at
##        k = 2 and k = 20: the difference lies within the equivalence margin
##        of +-.05, judged by its 95% Monte Carlo interval.
##    5.2 Part B: JZS (lmBF) is computed precisely enough: its reported
##        proportional error is below .05 in every dataset (recomputed with
##        10^5 and then 10^6 draws where it is above .02).
##    5.3 Part C: each prior's 95% intervals for x1 cover near nominally: in
##        every cell the coverage is above .93 or within 2 MCSE of .95.
##    5.4 Part C: each prior predicts at least as well as OLS: in every cell
##        its prediction-error ratio is below 1 + 2 MCSE.
##    5.5 Part C: the Stan model samples without divergent transitions.
##
##  Run from the project folder:  Rscript R/studies/05_behaviour.R
##  Smoke test (a corner of each design, 10 datasets per cell):
##    VALIDATION_SMOKE=1 Rscript R/studies/05_behaviour.R
##  Raw results are cached in cache/study_05/flatcor_<build>/ (delete it to
##  rerun); the tables and checks are saved to results/studies/study_05.rds.
##  Part C's fits of fit_regression() run 4 chains in parallel, so SimDesign
##  uses one worker per 4 cores there.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rsimsum)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # the methods, log_bayes_factors_one_coefficient()
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion()
source("R/data/data_sets.R")              # design_A(), generate_A(), with_stream_of()
source("R/data/data_set_A_factors.R")     # design_A_factors(), generate_A_factors()

study_folder <- study_cache_folder("05")
bootstrap_replicates <- if (smoke_mode) 100 else 1000

## Run a SimDesign simulation once and cache it, with its own resume file
cached_simulation <- function(name, ...) {
  file <- file.path(study_folder, paste0(name, ".rds"))
  if (file.exists(file)) { message("Using cached results: ", file); return(readRDS(file)) }
  simulation <- runSimulation(..., store_results = TRUE, progress = FALSE,
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = paste0("resume_", name, ".qs2")))
  saveRDS(simulation, file)
  simulation
}

## P(a value from higher exceeds one from lower), ties counting half: the
## Mann-Whitney form, from ranks
area_under_curve <- function(higher, lower) {
  ranks <- rank(c(higher, lower))
  (sum(ranks[seq_along(higher)]) - length(higher) * (length(higher) + 1) / 2) /
    (length(higher) * length(lower))
}


## =====================================================================
## PART A: ONE COEFFICIENT
## =====================================================================
design_one <- design_A()
if (smoke_mode) design_one <- subset(design_one, k %in% c(2, 20) & n == 50 & rho == 0 & others == "sparse")
n_sim_one <- n_sim_of("rate_per_cell")

## method names without spaces, for SimDesign's result columns
method_key <- setNames(make.names(model_methods), model_methods)

simulation_one <- cached_simulation("one_coefficient",
  design = design_one, replications = n_sim_one,
  generate = function(condition, fixed_objects) generate_A(condition, condition$REPLICATION),
  analyse = function(condition, dat, fixed_objects) {
    log_bfs <- log_bayes_factors_one_coefficient(dat$outcome, dat$predictors, tested_column = 1)
    setNames(log_bfs$log_BF10, method_key[log_bfs$method])   # log BF10 for x1, one per method
  },
  summarise = function(condition, results, fixed_objects)
    c(setNames(colMeans(results < -log(3)), paste0("null.", names(results))),
      setNames(colMeans(results > log(3)), paste0("effect.", names(results)))),
  parallel = number_of_cores > 1, ncores = number_of_cores, packages = c("flatcor", "BayesFactor"),
  control = list(include_replication_index = TRUE))

## every dataset's log BF10, one row per dataset and method
log_bf_one <- SimExtract(simulation_one, what = "results") |>
  pivot_longer(all_of(unname(method_key)), names_to = "key", values_to = "log_BF10") |>
  mutate(method = names(method_key)[match(key, method_key)],
         answer = case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                            TRUE ~ "undecided"))

## ---- Shares of each answer, with their Monte Carlo SE ---------------------------
shares_per_cell <- log_bf_one |>
  count(k, n, rho, others, slope_x1, method, answer, name = "datasets") |>
  complete(nesting(k, n, rho, others, slope_x1, method), answer, fill = list(datasets = 0)) |>
  group_by(k, n, rho, others, slope_x1, method) |>
  mutate(n_sim = sum(datasets), share = datasets / n_sim, mcse = mcse_of_proportion(share, n_sim)) |>
  ungroup()

## pooled over the cells at each k: the mean of the per-cell shares, whose
## MCSE follows from the per-cell MCSEs (independent cells)
shares_per_k <- shares_per_cell |>
  group_by(k, slope_x1, method, answer) |>
  summarise(cells = n(), share = mean(share), mcse = sqrt(sum(mcse^2)) / cells, .groups = "drop")

cat("\nShare of each answer at BF 3, pooled over the cells at each k (with MCSE):\n")
print(as.data.frame(shares_per_k |>
  mutate(share = sprintf("%.3f (%.3f)", share, mcse)) |> select(-mcse, -cells) |>
  pivot_wider(names_from = answer, values_from = share)), row.names = FALSE)


## ---- Discrimination: AUC null vs slope 0.1, and its difference from JZS ---------
## P(a value from higher exceeds one from lower), ties counting half: the
## Mann-Whitney form, from ranks, so the bootstrap stays fast at 625 per group

RNGkind("L'Ecuyer-CMRG"); set.seed(7)   # the bootstrap's own, fixed stream
auc <- log_bf_one |>
  filter(slope_x1 %in% c(0, 0.1)) |>
  group_by(k, n, rho, others) |>
  group_modify(function(cell, keys) {
    evidence_for_null <- function(method_name, slope)        # log BF01, in replication order
      -cell$log_BF10[cell$method == method_name & cell$slope_x1 == slope][
        order(cell$REPLICATION[cell$method == method_name & cell$slope_x1 == slope])]
    replications <- sum(cell$method == "JZS default" & cell$slope_x1 == 0)
    resamples <- replicate(bootstrap_replicates,
                           list(sample(replications, replace = TRUE), sample(replications, replace = TRUE)),
                           simplify = FALSE)
    auc_of <- function(method_name, rows_null = seq_len(replications), rows_effect = seq_len(replications))
      area_under_curve(evidence_for_null(method_name, 0)[rows_null],
                       evidence_for_null(method_name, 0.1)[rows_effect])
    bind_rows(lapply(model_methods, function(method_name) {
      versus_jzs <- vapply(resamples, function(rows)
        auc_of(method_name, rows[[1]], rows[[2]]) - auc_of("JZS default", rows[[1]], rows[[2]]), numeric(1))
      tibble(method = method_name, auc = auc_of(method_name),
             auc_minus_jzs = auc_of(method_name) - auc_of("JZS default"),
             mcse_minus_jzs = sd(versus_jzs),
             mcse_auc = sd(vapply(resamples, function(rows) auc_of(method_name, rows[[1]], rows[[2]]),
                                  numeric(1))))
    }))
  }) |>
  ungroup()

cat("\nAUC (null vs slope 0.1) minus JZS, mean over cells at each k:\n")
print(as.data.frame(auc |> group_by(k, method) |>
  summarise(auc = mean(auc), auc_minus_jzs = mean(auc_minus_jzs), .groups = "drop")), digits = 3)


## ---- Check 5.1: the share of true nulls confirmed at k = 2 and k = 20 -------------
null_confirmed <- function(method_name, k_value)
  filter(shares_per_k, method == method_name, k == k_value, slope_x1 == 0, answer == "says null")
k_drift <- bind_rows(lapply(model_methods, function(method_name) {
  at_2 <- null_confirmed(method_name, 2); at_20 <- null_confirmed(method_name, 20)
  difference <- at_20$share - at_2$share
  mcse <- sqrt(at_2$mcse^2 + at_20$mcse^2)
  tibble(method = method_name, at_k2 = at_2$share, at_k20 = at_20$share, difference = difference,
         lower = difference - 1.96 * mcse, upper = difference + 1.96 * mcse)
}))
for (row in seq_len(nrow(k_drift))) with(k_drift[row, ],
  record_check(sprintf("5.1 %s confirms true nulls at the same rate at k = 2 and k = 20", method),
               lower > -0.05 && upper < 0.05,
               sprintf("%.3f at k = 2, %.3f at k = 20: difference %.3f, 95%% Monte Carlo interval [%.3f, %.3f]; margin +-.05",
                       at_k2, at_k20, difference, lower, upper)))


## ---- Figure: what each test says, by k and true slope ---------------------------
outcome_shares <- shares_per_k |>
  filter(slope_x1 %in% c(0, 0.1, 0.2)) |>
  mutate(test = factor(method, levels = model_methods),
         answer = factor(answer, levels = c("says effect", "undecided", "says null")),
         truth = factor(paste("true slope", slope_x1), levels = paste("true slope", c(0, 0.1, 0.2))))
save_figure(
  ggplot(outcome_shares, aes(test, share, fill = answer)) +
    geom_col(width = 0.75) +
    geom_text(data = filter(outcome_shares, share >= 0.08), aes(label = sprintf("%.2f", share)),
              position = position_stack(vjust = 0.5), size = 2.7, colour = "white") +
    scale_fill_manual(values = c(`says null` = "#542788", undecided = "grey70", `says effect` = "#B35806"),
                      breaks = c("says null", "undecided", "says effect"), name = "at BF 3:") +
    facet_grid(truth ~ k, labeller = labeller(k = function(v) paste("k =", v))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    labs(x = NULL, y = "share of datasets",
         title = "What each Bayes factor says about x1 (the model with x1 against the model without it), as k grows",
         subtitle = sprintf("%d datasets per cell, pooled over n, correlation and pattern; correct: says null in the top row",
                            n_sim_one)),
  "study_05_one_coefficient_outcomes", width = 11, height = 8)


## =====================================================================
## PART B: A FACTOR
## =====================================================================
design_factor <- design_A_factors()
if (smoke_mode) design_factor <- subset(design_factor, q == 2 & n == 50 & other_predictors == 2)
r_squared_others <- 0.5        # what the other predictors explain (when the factor is null)

factor_methods <- c(names(contenders), "JZS (lmBF)", "Zellner-Siow")
factor_key <- setNames(make.names(factor_methods), factor_methods)

## The model with the factor against the model without it, every method; and
## the JZS proportional error
analyse_factor <- function(condition, dat, fixed_objects) {
  q <- condition$q; n <- condition$n; k <- q + condition$other_predictors
  others <- dat$others; group <- dat$group
  outcome <- outcome_A_factors(dat, condition, r_squared_others)
  design_full <- model.matrix(~ group + others)
  r_squared_full    <- sample_r_squared(design_full[, -1, drop = FALSE], outcome)
  r_squared_without <- sample_r_squared(others, outcome)
  joint <- function(method)
    log_bayes_factor_model(method, r_squared_full, n, k) -
      log_bayes_factor_model(method, r_squared_without, n, k - q)
  data_frame <- data.frame(y = outcome, group = group, others)
  formula_full    <- as.formula(paste("y ~ group +", paste(colnames(others), collapse = " + ")))
  formula_without <- as.formula(paste("y ~", paste(colnames(others), collapse = " + ")))
  ## lmBF integrates over two g's by Monte Carlo: give it the dataset's own
  ## stream, so its answer does not depend on the order or number of cores
  dataset_id <- paste(cell_id_A_factors(condition), condition$REPLICATION, sep = "-r")
  jzs <- with_stream_of(paste(dataset_id, "jzs", sep = "-"), {
    jzs_full <- BayesFactor::lmBF(formula_full, data_frame, progress = FALSE)
    for (iterations in c(1e5, 1e6))
      if (BayesFactor::extractBF(jzs_full)$error > 0.02)
        jzs_full <- BayesFactor::recompute(jzs_full, iterations = iterations, progress = FALSE)
    jzs_without <- BayesFactor::lmBF(formula_without, data_frame, progress = FALSE)
    c(log_bf = BayesFactor::extractBF(jzs_full, logbf = TRUE)$bf -
        BayesFactor::extractBF(jzs_without, logbf = TRUE)$bf,
      error = BayesFactor::extractBF(jzs_full)$error)
  })
  c(setNames(c(vapply(names(contenders), joint, numeric(1)), jzs[["log_bf"]], joint("Zellner-Siow")),
             unname(factor_key)),
    jzs_error = jzs[["error"]])
}

simulation_factor <- cached_simulation("factor",
  design = design_factor, replications = n_sim_one,
  generate = function(condition, fixed_objects) generate_A_factors(condition, condition$REPLICATION),
  analyse = analyse_factor,
  summarise = function(condition, results, fixed_objects) {
    log_bfs <- results[, unname(factor_key), drop = FALSE]
    c(setNames(colMeans(log_bfs < -log(3)), paste0("null.", names(log_bfs))),
      setNames(colMeans(log_bfs > log(3)), paste0("effect.", names(log_bfs))))
  },
  parallel = number_of_cores > 1, ncores = number_of_cores, packages = c("flatcor", "BayesFactor"),
  control = list(include_replication_index = TRUE))

raw_factor <- SimExtract(simulation_factor, what = "results")
jzs_errors <- raw_factor$jzs_error
log_bf_factor <- raw_factor |>
  select(-jzs_error) |>
  pivot_longer(all_of(unname(factor_key)), names_to = "key", values_to = "log_BF10") |>
  mutate(method = factor(names(factor_key)[match(key, factor_key)], levels = factor_methods),
         answer = case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                            TRUE ~ "undecided"))

shares_factor <- log_bf_factor |>
  count(q, n, other_predictors, factor_effect, method, answer, name = "datasets") |>
  complete(nesting(q, n, other_predictors, factor_effect, method), answer, fill = list(datasets = 0)) |>
  group_by(q, n, other_predictors, factor_effect, method) |>
  mutate(n_sim = sum(datasets), share = datasets / n_sim, mcse = mcse_of_proportion(share, n_sim)) |>
  ungroup()
shares_factor_over_q <- shares_factor |>
  group_by(n, other_predictors, factor_effect, method, answer) |>
  summarise(cells = n(), share = mean(share), mcse = sqrt(sum(mcse^2)) / cells, .groups = "drop")

cat("\nFactor: share of each answer at BF 3, averaged over q (with MCSE):\n")
print(as.data.frame(shares_factor_over_q |>
  mutate(share = sprintf("%.3f (%.3f)", share, mcse)) |> select(-mcse, -cells) |>
  pivot_wider(names_from = answer, values_from = share)), row.names = FALSE)

RNGkind("L'Ecuyer-CMRG"); set.seed(8)   # the bootstrap's own, fixed stream
auc_factor <- log_bf_factor |>
  filter(factor_effect %in% c("null", "small")) |>
  group_by(q, n, other_predictors) |>
  group_modify(function(cell, keys) {
    evidence_for_null <- function(method_name, effect)          # log BF01, in replication order
      -cell$log_BF10[cell$method == method_name & cell$factor_effect == effect][
        order(cell$REPLICATION[cell$method == method_name & cell$factor_effect == effect])]
    replications <- sum(cell$method == "JZS (lmBF)" & cell$factor_effect == "null")
    resamples <- replicate(bootstrap_replicates,
                           list(sample(replications, replace = TRUE), sample(replications, replace = TRUE)),
                           simplify = FALSE)
    auc_of <- function(method_name, rows_null = seq_len(replications), rows_small = seq_len(replications))
      area_under_curve(evidence_for_null(method_name, "null")[rows_null],
                       evidence_for_null(method_name, "small")[rows_small])
    bind_rows(lapply(factor_methods, function(method_name) {
      versus_jzs <- vapply(resamples, function(rows)
        auc_of(method_name, rows[[1]], rows[[2]]) - auc_of("JZS (lmBF)", rows[[1]], rows[[2]]), numeric(1))
      tibble(method = method_name, auc = auc_of(method_name),
             auc_minus_jzs = auc_of(method_name) - auc_of("JZS (lmBF)"), mcse_minus_jzs = sd(versus_jzs))
    }))
  }) |>
  ungroup()

cat("\nFactor: AUC (null vs small) minus JZS, mean over q, by n:\n")
print(as.data.frame(auc_factor |> group_by(method, n) |>
  summarise(auc_minus_jzs = mean(auc_minus_jzs), mcse = sqrt(sum(mcse_minus_jzs^2)) / n(),
            .groups = "drop")), digits = 3, row.names = FALSE)

record_check("5.2 JZS (lmBF) for a factor: proportional error < .05 in every dataset",
             max(jzs_errors) < 0.05,
             sprintf("largest %.3f over %d datasets", max(jzs_errors), length(jzs_errors)))


## =====================================================================
## PART C: ESTIMATION
## =====================================================================
design_estimation <- expand.grid(k = c(2, 10), n = c(50, 200), rho = c(0, 0.7),
                      others = c("sparse", "dense"), slope_x1 = 0.2, stringsAsFactors = FALSE)
if (smoke_mode) design_estimation <- subset(design_estimation, k == 2 & n == 50 & rho == 0)
n_sim_estimation <- n_sim_of("coverage_per_cell")
true_slope_x1 <- 0.2

## Stan settings: 4 chains of 1,000 draws; short chains in smoke mode
stan_chains   <- if (smoke_mode) 2 else 4
stan_warmup   <- if (smoke_mode) 200 else 1000
stan_sampling <- if (smoke_mode) 200 else 1000
bayesfactor_iterations <- if (smoke_mode) 1000 else 4000
workers <- max(1, number_of_cores %/% stan_chains)

estimation_methods <- c(names(contenders), "JZS default", "Zellner-Siow", "ordinary (lm)")
estimation_key <- setNames(make.names(estimation_methods), estimation_methods)


## ---- The simulation ------------------------------------------------------------
simulate_data <- function(condition, fixed_objects) generate_A(condition, condition$REPLICATION)

analyse_dataset <- function(condition, dat, fixed_objects) {
  dataset_id <- paste(cell_id_A(condition), condition$REPLICATION, sep = "-r")
  k <- condition$k
  covariance  <- matrix(condition$rho, k, k); diag(covariance) <- 1
  true_slopes <- slopes_A(k, condition$others, condition$slope_x1)
  data <- data.frame(dat$predictors, y = dat$outcome)
  formula <- reformulate(colnames(dat$predictors), response = "y")
  ## every fit's seed follows from its dataset's ID, so results do not depend
  ## on the number of cores or the order of the fits
  seed_of <- function(label) abs(digest::digest2int(paste(dataset_id, label))) %% .Machine$integer.max

  summary_of <- function(draws) {      # posterior mean and 95% interval of x1, and the excess prediction error
    estimates <- colMeans(draws)
    error <- estimates - true_slopes
    c(estimate = estimates[["x1"]], sd = sd(draws[, "x1"]),
      lower = quantile(draws[, "x1"], 0.025, names = FALSE),
      upper = quantile(draws[, "x1"], 0.975, names = FALSE),
      prediction_error = drop(t(error) %*% covariance %*% error))
  }
  per_method <- list()
  divergent  <- 0
  for (method in names(contenders)) {
    fit <- suppressWarnings(fit_regression(formula, data, prior = contenders[[method]],
                                           chains = stan_chains, iter_warmup = stan_warmup,
                                           iter_sampling = stan_sampling, seed = seed_of(method)))
    divergent <- divergent + fit$divergent
    per_method[[method]] <- summary_of(as.matrix(as.data.frame(fit$draws)[, colnames(dat$predictors),
                                                                         drop = FALSE]))
  }
  bayesfactor_draws <- function(scale, label)
    with_stream_of(paste(dataset_id, label), as.matrix(as.data.frame(
      BayesFactor::lmBF(formula, data, rscaleCont = scale, posterior = TRUE,
                        iterations = bayesfactor_iterations, progress = FALSE))[, colnames(dat$predictors),
                                                                               drop = FALSE]))
  per_method[["JZS default"]]  <- summary_of(bayesfactor_draws(sqrt(2) / 4, "JZS"))
  per_method[["Zellner-Siow"]] <- summary_of(bayesfactor_draws(1, "ZS"))
  ols_fit <- lm(formula, data)
  ols_error <- coef(ols_fit)[colnames(dat$predictors)] - true_slopes
  per_method[["ordinary (lm)"]] <- c(estimate = coef(ols_fit)[["x1"]],
                                     sd = summary(ols_fit)$coefficients["x1", "Std. Error"],
                                     lower = confint(ols_fit)["x1", 1], upper = confint(ols_fit)["x1", 2],
                                     prediction_error = drop(t(ols_error) %*% covariance %*% ols_error))
  values <- unlist(lapply(estimation_methods, function(method)
    setNames(per_method[[method]], paste(names(per_method[[method]]), estimation_key[[method]], sep = "."))))
  c(values, divergent = divergent)
}

## per cell: the divergent transitions of the Stan fits (the estimands are
## computed from the stored results below)
summarise_cell <- function(condition, results, fixed_objects) c(divergent = sum(results[, "divergent"]))

simulation_estimation <- cached_simulation("estimation",
  design = design_estimation, replications = n_sim_estimation, generate = simulate_data,
  analyse = analyse_dataset, summarise = summarise_cell,
  parallel = workers > 1, ncores = workers, packages = c("flatcor", "BayesFactor"),
  control = list(include_replication_index = TRUE))

## one row per dataset and method: x1's estimate and interval, and the excess
## prediction error
per_dataset <- SimExtract(simulation_estimation, what = "results")
estimates <- per_dataset |>
  select(-divergent) |>
  pivot_longer(matches("^(estimate|sd|lower|upper|prediction_error)[.]"),
               names_to = c(".value", "key"), names_pattern = "^([a-z_]+)[.](.*)$") |>
  mutate(method = names(estimation_key)[match(key, estimation_key)],
         cell = sprintf("k = %d, n = %d, rho = %.1f, %s", k, n, rho, others))


## ---- Bias and coverage of x1, with their Monte Carlo SE (rsimsum) -------------
x1_summary <- simsum(as.data.frame(rename(estimates, ci_lower = lower, ci_upper = upper)),
                     estvarname = "estimate", se = "sd", true = true_slope_x1, methodvar = "method",
                     by = c("k", "n", "rho", "others"), ci.limits = c("ci_lower", "ci_upper"))
## rbias: the bias as a share of the true value 0.2 (rsimsum's relative bias)
bias_and_coverage <- tidy(summary(x1_summary, stats = c("rbias", "cover"))) |>
  mutate(stat = recode(stat, rbias = "bias"),
         across(c(k, n, rho), ~ as.numeric(as.character(.x))))   # rsimsum returns the by-variables as factors

cat("\nBias of x1 (share of its true value 0.2) and coverage of its 95% interval, with MCSE:\n")
print(as.data.frame(bias_and_coverage |>
  mutate(value = sprintf("%.3f (%.3f)", est, mcse)) |>
  select(stat, k, n, rho, others, method, value) |>
  pivot_wider(names_from = method, values_from = value)), row.names = FALSE)


## ---- Excess prediction error relative to OLS, with a paired bootstrap SE --------
RNGkind("L'Ecuyer-CMRG"); set.seed(9)   # the bootstrap's own, fixed stream
prediction_ratio <- estimates |>
  group_by(k, n, rho, others) |>
  group_modify(function(cell, keys) {
    error_of <- function(method_name)                       # in replication order
      cell$prediction_error[cell$method == method_name][order(cell$REPLICATION[cell$method == method_name])]
    replications <- sum(cell$method == "ordinary (lm)")
    resamples <- replicate(bootstrap_replicates, sample(replications, replace = TRUE), simplify = FALSE)
    ols <- error_of("ordinary (lm)")
    bind_rows(lapply(setdiff(estimation_methods, "ordinary (lm)"), function(method_name) {
      errors <- error_of(method_name)
      tibble(method = method_name, ratio = mean(errors) / mean(ols),
             mcse = sd(vapply(resamples, function(rows) mean(errors[rows]) / mean(ols[rows]), numeric(1))))
    }))
  }) |>
  ungroup()

cat("\nExcess prediction error relative to OLS (below 1: better than OLS), with MCSE:\n")
print(as.data.frame(prediction_ratio |>
  mutate(value = sprintf("%.3f (%.3f)", ratio, mcse)) |> select(-ratio, -mcse) |>
  pivot_wider(names_from = method, values_from = value)), row.names = FALSE)


## ---- Checks 5.3 to 5.5 -----------------------------------------------------------
for (prior_name in names(contenders)) {
  coverage <- filter(bias_and_coverage, stat == "cover", method == prior_name)
  record_check(sprintf("5.3 %s: 95%% intervals for x1 cover near nominally", prior_name),
               all(coverage$est > 0.93 | abs(coverage$est - 0.95) <= 2 * coverage$mcse),
               sprintf("lowest coverage %.3f (MCSE %.3f) over %d cells",
                       min(coverage$est), coverage$mcse[which.min(coverage$est)], nrow(coverage)))
}
for (prior_name in names(contenders)) {
  ratio <- filter(prediction_ratio, method == prior_name)
  record_check(sprintf("5.4 %s predicts at least as well as OLS", prior_name),
               all(ratio$ratio < 1 + 2 * ratio$mcse),
               sprintf("highest ratio to OLS %.3f (MCSE %.3f) over %d cells",
                       max(ratio$ratio), ratio$mcse[which.max(ratio$ratio)], nrow(ratio)))
}
record_check("5.5 the Stan model samples without divergent transitions",
             sum(per_dataset$divergent) == 0,
             sprintf("%d divergent transitions over %d fits (%d datasets x 3 priors)",
                     sum(per_dataset$divergent), 3 * nrow(per_dataset), nrow(per_dataset)))


## ---- Figures -------------------------------------------------------------------
cells_facet <- facet_grid(others + rho ~ k,
  labeller = labeller(k = function(v) paste("k =", v), rho = function(v) paste("correlation =", v)))
coverage_plot <- filter(bias_and_coverage, stat == "cover")
save_figure(
  ggplot(coverage_plot, aes(n, est, colour = method)) +
    geom_hline(yintercept = 0.95, linetype = "dotted") +
    geom_pointrange(aes(ymin = est - 2 * mcse, ymax = est + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "coverage of x1's 95% interval",
         title = "Do the intervals contain the true slope of x1 (0.2)?",
         subtitle = sprintf("%d datasets per cell; bars: 2 Monte Carlo SE", n_sim_estimation)),
  "study_05_estimation_coverage", width = 9, height = 10)

bias_plot <- filter(bias_and_coverage, stat == "bias")
save_figure(
  ggplot(bias_plot, aes(n, est, colour = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_pointrange(aes(ymin = est - 2 * mcse, ymax = est + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "bias of x1, as a share of its true value",
         title = "Bias: how far estimates of x1 are pulled towards zero",
         subtitle = sprintf("%d datasets per cell; bars: 2 Monte Carlo SE", n_sim_estimation)),
  "study_05_estimation_bias", width = 9, height = 10)

save_figure(
  ggplot(prediction_ratio, aes(n, ratio, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_pointrange(aes(ymin = ratio - 2 * mcse, ymax = ratio + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) +
    method_scales(setdiff(estimation_methods, "ordinary (lm)"), linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "excess prediction error / OLS",
         title = "Does the shrinkage pay off for prediction?",
         subtitle = "below 1: predicts new observations better than ordinary regression; bars: 2 Monte Carlo SE"),
  "study_05_estimation_prediction_error", width = 9, height = 10)

save_study_results("05", list(shares_per_cell = shares_per_cell, shares_per_k = shares_per_k, auc = auc,
                              k_drift = k_drift, factor_shares = shares_factor,
                              factor_shares_over_q = shares_factor_over_q, factor_auc = auc_factor,
                              bias_and_coverage = bias_and_coverage, prediction_ratio = prediction_ratio,
                              divergent = sum(per_dataset$divergent)))
report_checks()
