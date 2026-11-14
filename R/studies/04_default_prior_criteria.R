## =====================================================================
##  Study 4: does each prior meet the default-prior criteria?
##  (docs/appendix/C_study_design.md, study 4)
##
##  Aims. The criteria of Bayarri, Berger, Forte and Garcia-Donato (2012) for
##  each of the three priors: predictive matching, information consistency,
##  intrinsic consistency, model-selection consistency for a fixed and a
##  growing number of predictors k, and invariance under rescaling and
##  recoding of the predictors.
##  Data. None for the formula checks; simulated R-squared values and datasets
##  for the rates (per-draw streams from their IDs); data set D for growing k
##  (R/data/data_set_D.R).
##  Estimands.
##    Fixed k: the slope of the median log Bayes factor for the truth against
##    log n. Expected (asymptotically): a for the intercept-only model, with a
##    the shape of R-squared ~ Beta(a, 1) (1/2, 1 or q/2); q/2 when the
##    smaller model has a real effect and the larger adds q null predictors.
##    Growing k: the same slope, which must be positive (consistent).
##  Methods. The three priors (whole-model Bayes factors); JZS and Zellner-Siow
##  as references for growing k.
##  Performance. Slopes and medians with bootstrap MCSE (resampling the draws
##  within each n). The expected slopes are asymptotic, so a slope passes when
##  its 95% Monte Carlo interval lies within the stated tolerance of the
##  expected value (the tolerance covers finite-n deviation; the interval
##  covers Monte Carlo error).
##  n_sim. 2,000 R-squared draws per setting for the rates and for growing k
##  under the null (n_sim_of("rate_draws"); slope MCSE about 0.02); 400 datasets per
##  setting for growing k with a true model (n_sim_of("growing_k_data")).
##
##  Checks (each for each prior):
##    4.1 predictive matching: BF = 1 at n = k + 1
##    4.2 information consistency: evidence increases without bound as
##        R-squared -> 1
##    4.3 information consistency starts exactly at n = k + 3
##    4.4 consistency under the null: one null predictor next to two real
##        ones, 95% interval of the slope within (.35, .65)
##    4.5 the rate with a real effect in the smaller model: 95% interval of
##        the slope within 20% of q/2, q = 1, 3, 6
##    4.6 the rate against the intercept-only model: 95% interval within 30%
##        of a, q = 1, 3, 6
##    4.7 consistency under the alternative: median log BF10 for a real effect
##        increases with n, and its 95% lower bound at n = 5,000 exceeds 20
##    4.8 intrinsic consistency: the prior with the sample covariance is
##        indistinguishable from the prior with the population covariance at
##        n = 100,000 and distinguishable at n = 30 (two-sample KS against its
##        .01 critical value)
##    4.9 growing k, null true: evidence for the null grows; 95% lower bound
##        of the slope above 0.1
##    4.10 growing k, small model true: evidence for it grows; 95% lower
##        bound of the slope above 0.1
##    4.11 invariance: rescaling, shifting and recoding the predictors (and
##        the outcome) leaves the Bayes factor unchanged, within a relative
##        1e-8
##
##  Run from the project folder:  Rscript R/studies/04_default_prior_criteria.R
##  Smoke test:  VALIDATION_SMOKE=1 Rscript R/studies/04_default_prior_criteria.R
##  Raw results are cached in cache/study_04/flatcor_<build>/ (delete it to rerun); the
##  tables and checks are saved to results/studies/study_04.rds.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # contenders, model_methods, log_bayes_factor_model()
source("R/shared/precision.R")            # n_sim_of(), smoke_mode
source("R/data/data_sets.R")            # with_stream_of()
source("R/data/data_set_D.R")           # design_D(), generate_D_null(), generate_D_true()

study_folder <- study_cache_folder("04")

priors <- c("flat_r", "flat_partial_r", "flat_r_squared")
log_bf <- function(prior, r_squared, n, k)
  if (k == 0) 0 else log_bayes_factor_vs_intercept_only(r_squared, n, k, prior)
bootstrap_replicates <- if (smoke_mode) 100 else 1000

## Run a SimDesign simulation once and cache it
cached_simulation <- function(name, ...) {
  file <- file.path(study_folder, paste0(name, ".rds"))
  if (file.exists(file)) { message("Using cached results: ", file); return(readRDS(file)) }
  ## SimDesign's resume file goes to this study's cache folder, one per
  ## simulation: its default (named after the machine, in the working
  ## directory) collides with other simulations running at the same time
  simulation <- runSimulation(..., parallel = number_of_cores > 1, ncores = number_of_cores,
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = paste0("resume_", name, ".qs2")),
                              packages = c("flatcor", "BayesFactor"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE)
  saveRDS(simulation, file)
  simulation
}

## Slope of the median of value against log n, with a bootstrap MCSE from
## resampling the draws within each n (a data frame with columns n, value)
slope_of_medians <- function(draws) {
  by_n <- split(draws$value, draws$n)
  sizes <- as.numeric(names(by_n))
  slope <- function(medians) coef(lm(medians ~ log(sizes)))[[2]]
  estimate <- slope(vapply(by_n, median, numeric(1)))
  replicates <- replicate(bootstrap_replicates,
                          slope(vapply(by_n, function(v) median(sample(v, replace = TRUE)), numeric(1))))
  mcse <- sd(replicates)
  c(slope = estimate, mcse = mcse, lower = estimate - 1.96 * mcse, upper = estimate + 1.96 * mcse)
}


## =====================================================================
##  Part 1: the formula checks (no simulation)
## =====================================================================
for (prior in priors) {
  ## 4.1 predictive matching. At n = k + 1 the package returns BF = 1 without
  ## integrating (every model fits exactly), so the check computes the integral
  ## itself: the package's likelihood ratio averaged over R-squared ~ Beta(a, 1)
  ## by integrate(), for a sample R-squared just below 1.
  ks <- 1:6
  matching <- sapply(ks, function(k) {
    n <- k + 1; a <- prior_shape(prior, k)
    integrand <- function(f) exp(log_likelihood_ratio_block(f, 1 - 1e-12, n, k, k)) *
      2 * a * f^(2 * a - 1) / (1 + f^2)^(a + 1)
    sum(sapply(list(c(0, 1), c(1, 100), c(100, Inf)), function(range)
      integrate(integrand, range[1], range[2], rel.tol = 1e-10)$value))
  })
  record_check(sprintf("4.1 %s: predictive matching, BF = 1 at n = k + 1 (k = %d-%d)", prior, min(ks), max(ks)),
               all(abs(matching - 1) < 1e-3) && all(exp(sapply(ks, function(k) log_bf(prior, 1 - 1e-12, k + 1, k))) == 1),
               paste(sprintf("%.4f", matching), collapse = " "))

  ## 4.2 information consistency: strictly increasing, above 50 at R2 = .9999
  information <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(r2) log_bf(prior, r2, 20, 3))
  record_check(sprintf("4.2 %s: information consistency", prior),
               all(diff(information) > 0) && information[4] > 50,
               paste("log BF at R2 = .9 ... .9999:", paste(sprintf("%.1f", information), collapse = " ")))

  ## 4.3 the exact boundary: bounded at n = k + 2, growing at n = k + 3
  boundary <- sapply(1:3, function(k) {
    at <- function(n) sapply(c(1 - 1e-6, 1 - 1e-8), function(r2) log_bf(prior, r2, n, k))
    c(bounded_at_k_plus_2 = diff(at(k + 2)) < 0.01, growing_at_k_plus_3 = diff(at(k + 3)) > 0.2)
  })
  record_check(sprintf("4.3 %s: information consistency starts exactly at n = k + 3 (k = 1-3)", prior),
               all(boundary), "bounded at n = k + 2, growing at n = k + 3")
}


## =====================================================================
##  Part 2: the rates for a fixed number of predictors (simulated)
## =====================================================================
## Each setting compares a smaller (true) model with a larger one. The extra
## columns of the larger model are null, so their partial R-squared given the
## smaller model is exactly Beta(q/2, (n - k_small - 1 - q)/2) and is drawn
## directly; the smaller model's R-squared comes from simulated data when it
## has real effects. "alternative" is the reverse: a real x1 = 0.15 next to
## x2 = 0.2 and a null x3, evidence for the model with x1.
rates_settings <- bind_rows(
  tibble(setting = "x3 null next to two real",  k_small = 2, q = 1),
  tibble(setting = "intercept only",            k_small = 0, q = c(1, 3, 6)),
  tibble(setting = "x1 real",                   k_small = 1, q = c(1, 3, 6)),
  tibble(setting = "alternative",               k_small = 2, q = 1))
rates_design <- rates_settings |>
  crossing(n = c(50, 200, 1000, 5000)) |>
  filter(setting != "alternative" | n >= 200)
if (smoke_mode) rates_design <- filter(rates_design, q == 1, n %in% c(200, 1000))
rates_design <- as.data.frame(rates_design)

simulate_rates <- function(condition, fixed_objects) {
  n <- condition$n; q <- condition$q
  with_stream_of(sprintf("R5-%s-q%d-n%d-r%d", condition$setting, q, n, condition$REPLICATION), {
    r_squared_of <- function(x, y) sample_r_squared(x, y)
    if (condition$setting == "alternative") {
      x <- matrix(rnorm(n * 3), n); y <- drop(x %*% c(0.15, 0.2, 0)) + rnorm(n)
      return(c(r2_small = r_squared_of(x[, 2:3], y), r2_large = r_squared_of(x, y)))
    }
    r2_small <- switch(condition$setting,
      `intercept only` = 0,
      `x1 real` = { x <- rnorm(n); r_squared_of(matrix(x), 0.5 * x + rnorm(n)) },
      `x3 null next to two real` = { x <- matrix(rnorm(n * 2), n)
                                     r_squared_of(x, drop(x %*% c(0.3, 0.2)) + rnorm(n)) })
    partial <- rbeta(1, q / 2, (n - condition$k_small - 1 - q) / 2)
    c(r2_small = r2_small, r2_large = r2_small + (1 - r2_small) * partial)
  })
}

## log BF for the truth: the smaller model against the larger, except for
## "alternative", where the larger model (with the real x1) is the truth
analyse_rates <- function(condition, dat, fixed_objects) {
  k_small <- condition$k_small; k_large <- k_small + condition$q
  sapply(setNames(priors, priors), function(prior) {
    small_vs_large <- log_bf(prior, dat[["r2_small"]], condition$n, k_small) -
                      log_bf(prior, dat[["r2_large"]], condition$n, k_large)
    if (condition$setting == "alternative") -small_vs_large else small_vs_large
  })
}
summarise_medians <- function(condition, results, fixed_objects) apply(results, 2, median)

rates <- cached_simulation("rates", design = rates_design, replications = n_sim_of("rate_draws"),
                           generate = simulate_rates, analyse = analyse_rates,
                           summarise = summarise_medians)
rate_draws <- SimExtract(rates, what = "results") |>
  pivot_longer(all_of(priors), names_to = "prior", values_to = "value")

RNGkind("L'Ecuyer-CMRG"); set.seed(5)   # the bootstrap's own, fixed stream
rate_slopes <- rate_draws |>
  filter(setting != "alternative") |>
  group_by(setting, q, prior) |>
  group_modify(function(d, keys) as_tibble(as.list(slope_of_medians(d)))) |>
  ungroup() |>
  mutate(expected = case_when(setting == "intercept only" ~ mapply(prior_shape, prior, q),
                              TRUE ~ q / 2))
cat("\nSlope of the median log BF for the truth against log n (with MCSE):\n")
print(as.data.frame(rate_slopes), digits = 3, row.names = FALSE)

interval_within <- function(rows, tolerance)
  all(rows$lower > rows$expected * (1 - tolerance) & rows$upper < rows$expected * (1 + tolerance))
describe <- function(rows)
  paste(sprintf("q = %d: %.2f [%.2f, %.2f] (expected %.2f)", rows$q, rows$slope, rows$lower, rows$upper,
                rows$expected), collapse = "; ")

for (p in priors) {
  ## 4.4 one null predictor next to two real ones: the slope is near 1/2
  rows <- filter(rate_slopes, prior == p, setting == "x3 null next to two real")
  record_check(sprintf("4.4 %s: consistency under the null, slope within (.35, .65)", p),
               nrow(rows) == 0 || (rows$lower > 0.35 && rows$upper < 0.65), describe(rows))
  ## 4.5 a real effect in the smaller model: the slope is q/2
  rows <- filter(rate_slopes, prior == p, setting == "x1 real")
  record_check(sprintf("4.5 %s: smaller model x1, slope within 20%% of q/2", p),
               interval_within(rows, 0.2), describe(rows))
  ## 4.6 the intercept-only model: the slope is a, the shape of the prior
  rows <- filter(rate_slopes, prior == p, setting == "intercept only")
  record_check(sprintf("4.6 %s: smaller model intercept only, slope within 30%% of a", p),
               interval_within(rows, 0.3), describe(rows))
}

## 4.7 a real effect: the median log BF10 grows with n, above 20 at n = 5,000
RNGkind("L'Ecuyer-CMRG"); set.seed(57)   # the bootstrap's own, fixed stream
alternative_medians <- rate_draws |>
  filter(setting == "alternative") |>
  group_by(prior, n) |>
  summarise(median = median(value),
            mcse = sd(replicate(bootstrap_replicates, median(sample(value, replace = TRUE)))),
            .groups = "drop")
for (p in priors) {
  rows <- arrange(filter(alternative_medians, prior == p), n)
  largest <- tail(rows, 1)
  record_check(sprintf("4.7 %s: consistency under the alternative", p),
               nrow(rows) >= 2 && all(diff(rows$median) - 1.96 * sqrt(head(rows$mcse, -1)^2 +
                                                                       tail(rows$mcse, -1)^2) > 0) &&
                 (smoke_mode || (largest$n == 5000 && largest$median - 1.96 * largest$mcse > 20)),
               paste(sprintf("n = %d: %.1f (MCSE %.2f)", rows$n, rows$median, rows$mcse), collapse = "; "))
}


## =====================================================================
##  Part 3: intrinsic consistency
## =====================================================================
## The prior on the coefficients uses the sample covariance of the
## predictors, so it depends on n; it must tend to the prior with the
## population covariance. Three predictors with correlation .5; the implied
## prior of the first coefficient (noise units) with the sample covariance
## against the one with the population covariance, by a two-sample KS
## distance against its .01 critical value 1.63 sqrt(2 / draws).
intrinsic_draws <- n_sim_of("ks_draws")
critical_two_sample <- 1.63 * sqrt(2 / intrinsic_draws)
population <- matrix(0.5, 3, 3); diag(population) <- 1
for (p in priors) {
  limit_draws <- with_stream_of(paste0("I5-", p, "-limit"),
                                draw_coefficients_from_prior(intrinsic_draws, population, prior = p)[, 1])
  distance <- sapply(c(30, 300, 100000), function(n) with_stream_of(sprintf("I5-%s-n%d", p, n), {
    x <- matrix(rnorm(n * 3), n) %*% chol(population)
    sample_covariance <- crossprod(scale(x, scale = FALSE)) / n
    suppressWarnings(unname(ks.test(draw_coefficients_from_prior(intrinsic_draws, sample_covariance,
                                                                 prior = p)[, 1], limit_draws)$statistic))
  }))
  record_check(sprintf("4.8 %s: intrinsic consistency (prior with Sigma_n -> prior with Sigma)", p),
               distance[3] < critical_two_sample && distance[1] > critical_two_sample,
               sprintf("KS distance at n = 30, 300, 100,000: %s (critical value %.4f, %d draws)",
                       paste(sprintf("%.4f", distance), collapse = " "), critical_two_sample,
                       intrinsic_draws))
}


## =====================================================================
##  Part 4: k growing with n (data set D)
## =====================================================================
growing_design <- design_D()
if (smoke_mode) growing_design <- subset(growing_design, n %in% c(50, 100) & growth == "k = n / 4")
growing_design <- as.data.frame(growing_design)

## null true: log BF01 of the null against the k-predictor model
null_draws <- cached_simulation("growing_null", design = growing_design,
  replications = n_sim_of("rate_draws"),
  generate = function(condition, fixed_objects) generate_D_null(condition, condition$REPLICATION),
  analyse = function(condition, dat, fixed_objects) {
    k <- k_of_D(condition$n, condition$growth)
    sapply(setNames(model_methods, make.names(model_methods)), function(method)
      -log_bayes_factor_model(method, dat, condition$n, k))
  },
  summarise = summarise_medians)

## small model true: log BF of the true 3-predictor model against the full one
true_draws <- cached_simulation("growing_true", design = growing_design,
  replications = n_sim_of("growing_k_data"),
  generate = function(condition, fixed_objects) generate_D_true(condition, condition$REPLICATION),
  analyse = function(condition, dat, fixed_objects) {
    k <- k_of_D(condition$n, condition$growth)
    r2_true <- sample_r_squared(dat$predictors[, 1:3], dat$outcome)
    r2_full <- sample_r_squared(dat$predictors, dat$outcome)
    sapply(setNames(model_methods, make.names(model_methods)), function(method)
      log_bayes_factor_model(method, r2_true, condition$n, 3) -
        log_bayes_factor_model(method, r2_full, condition$n, k))
  },
  summarise = summarise_medians)

method_of_key <- setNames(model_methods, make.names(model_methods))
growing <- bind_rows(
  SimExtract(null_draws, what = "results") |> mutate(question = "null true: log BF(null vs k-predictor model)"),
  SimExtract(true_draws, what = "results") |> mutate(question = "3-predictor model true: log BF(true vs full)")) |>
  pivot_longer(all_of(names(method_of_key)), names_to = "key", values_to = "value") |>
  mutate(method = method_of_key[key])

RNGkind("L'Ecuyer-CMRG"); set.seed(59)   # the bootstrap's own, fixed stream
growing_slopes <- growing |>
  group_by(question, growth, method) |>
  group_modify(function(d, keys) as_tibble(as.list(slope_of_medians(d)))) |>
  ungroup()
growing_medians <- growing |>
  group_by(question, growth, n, method) |>
  summarise(median = median(value), .groups = "drop")
cat("\nGrowing k: slope of the median log BF for the truth against log n (with MCSE):\n")
print(as.data.frame(growing_slopes), digits = 3, row.names = FALSE)

for (p in names(contenders)) for (g in unique(growing$growth)) {
  for (check in list(list(id = "4.9", question = "null true: log BF(null vs k-predictor model)",
                          label = "null true: evidence for the null grows"),
                     list(id = "4.10", question = "3-predictor model true: log BF(true vs full)",
                          label = "small model true: evidence for it grows"))) {
    row <- filter(growing_slopes, method == p, growth == g, question == check$question)
    medians <- arrange(filter(growing_medians, method == p, growth == g, question == check$question), n)
    record_check(sprintf("%s %s, %s, %s", check$id, p, g, check$label),
                 row$lower > 0.1,
                 sprintf("slope %.2f [%.2f, %.2f]; median log BF at n = %s: %s", row$slope, row$lower,
                         row$upper, paste(medians$n, collapse = ", "),
                         paste(sprintf("%.1f", medians$median), collapse = " ")))
  }
}

save_figure(
  ggplot(growing_medians, aes(n, median, colour = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = sort(unique(growing_medians$n))) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 3, 10, 30, 100, 300)) +
    method_scales(model_methods) +
    facet_grid(question ~ growth) +
    labs(x = "n (log scale)", y = "median log Bayes factor for the truth (pseudo-log scale)",
         title = "Consistency when the number of predictors grows with n",
         subtitle = sprintf("%d R-squared draws (null) and %d datasets (true model) per setting",
                            n_sim_of("rate_draws"), n_sim_of("growing_k_data"))),
  "study_04_growing_k", width = 11, height = 7)


## ---- 4.11 invariance ----------------------------------------------------------------
## The Bayes factor of each predictor (the model with it against the model
## without it) and of the whole model, before and after rescaling, shifting
## and recoding the predictors and the outcome
invariance_data <- local({
  set.seed(42)
  correlations <- matrix(0.4, 3, 3); diag(correlations) <- 1
  x <- matrix(rnorm(60 * 3), 60) %*% chol(correlations)
  data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], y = drop(x %*% c(0.3, 0.2, 0)) + rnorm(60))
})
bayes_factors_of <- function(data, formula_full, prior) {
  full <- lm(formula_full, data)
  terms_full <- attr(terms(full), "term.labels")
  c(whole = bayes_factor_models(full, lm(y ~ 1, data), prior = prior),
    sapply(terms_full, function(term)
      bayes_factor_models(full, update(full, paste(". ~ . -", term)), prior = prior)))
}
for (p in priors) {
  original <- bayes_factors_of(invariance_data, y ~ x1 + x2 + x3, p)
  rescaled <- bayes_factors_of(transform(invariance_data, x1 = 100 * x1 + 7, x2 = -0.01 * x2,
                                         x3 = x3 - 50, y = 1000 * y - 3), y ~ x1 + x2 + x3, p)
  recoded_whole <- bayes_factors_of(transform(invariance_data, sum_12 = x1 + x2, difference_12 = x1 - x2),
                                    y ~ sum_12 + difference_12 + x3, p)[["whole"]]
  largest <- max(abs(rescaled / original - 1), abs(recoded_whole / original[["whole"]] - 1))
  record_check(sprintf("4.11 %s: invariance under rescaling, shifting and recoding", p), largest < 1e-8,
               sprintf("largest relative difference %.1e", largest))
}

save_study_results("04", list(rate_slopes = rate_slopes, alternative_medians = alternative_medians,
                              growing_slopes = growing_slopes, growing_medians = growing_medians))
report_checks()
