## =====================================================================
##  Study 2: are the Bayes factors computed correctly?
##  (docs/appendix/C_study_design.md, study 2)
##
##  Aims. flatcor computes every Bayes factor between two models from one
##  route, the same for all three priors: the whole-model Bayes factor against
##  the intercept-only model, the likelihood ratio of the model's columns
##  averaged over R-squared ~ Beta(a, 1) by integrate() (flat_r a = 1/2,
##  flat_partial_r a = k/2, flat_r_squared a = 1). Does it equal the integral
##  it stands for, for every prior?
##  Data. Fixed grids of summary statistics (R-squared, n, k), extremes
##  included, and three simulated datasets for the routes that need data.
##  Estimands. The log Bayes factor against the intercept-only model; the log
##  likelihood ratio inside it.
##  Methods. The package, against routes that share no code with it, the same
##  four for every prior:
##    2.1 adaptive integration (stats::integrate()) of a second derivation:
##        the prior as a mixture of g-priors (for flat_r, averaged over
##        directions; for flat_r_squared with one predictor, which is not a
##        g-prior mixture, over the effect itself);
##    2.2 the likelihood ratio the three priors share, against its exact
##        hypergeometric series summed in 200-bit arithmetic (Rmpfr);
##    2.3 brute force: draw coefficients from the prior, average how well they
##        predict (no formula);
##    2.4 bridge sampling (bridgesampling) on draws from the package's Stan
##        model.
##  Performance. The absolute log difference from the reference. For the
##  Monte Carlo routes (2.3, 2.4) it is judged against their own Monte Carlo
##  error, and each runs until that error is below 0.01 on the log scale (a
##  stated precision for an estimator, not a stopping rule on the result).
##  n_sim. Deterministic, except the Monte Carlo routes (to MCSE 0.01).
##
##  Checks, for each prior unless stated:
##    2.1a integrate() of the g-prior form, ordinary evidence (|log BF| < 50):
##         within 1e-5
##    2.1b the same, all settings: within a relative error of 1e-6 (|log BF|
##         above 1) or 1e-5
##    2.2  the likelihood ratio (prior-free, shared by the three priors) against
##         its 200-bit series, q = 1 and q = k: within 1e-6
##    2.3  brute force, run to MCSE 0.01: within 3 MCSE
##    2.4  bridge sampling on the package's Stan fit, run to an estimated error
##         of 0.01: within 3 of its errors. bridgesampling's error allows for
##         autocorrelation (Fruhwirth-Schnatter, 2004), not for a biased
##         sampler, so the fits use adapt_delta = 0.99.
##  Invariance under rescaling and recoding is a default-prior criterion, in
##  study 4.
##
##  Run from the project folder:  Rscript R/studies/02_bayes_factors.R
##  Smoke test (small grids, few draws, a few minutes):
##    VALIDATION_SMOKE=1 Rscript R/studies/02_bayes_factors.R
##  Needs cmdstanr with CmdStan (check 2.4), Rmpfr and bridgesampling.
##  Raw results are cached in cache/study_02/flatcor_<build>/ (delete it to
##  rerun); the checks are saved to results/studies/study_02.rds.
## =====================================================================
library(dplyr)
source("R/shared/settings.R")             # flatcor, cores, record_check(), parallel_map()
source("R/shared/precision.R")            # precision_targets, smoke_mode

study_folder <- study_cache_folder("02")

priors <- c("flat_r", "flat_partial_r", "flat_r_squared")
log_bf <- function(prior, r_squared, n, k) log_bayes_factor_vs_intercept_only(r_squared, n, k, prior)

## a cached computation: the expression is evaluated once and its value saved
cached <- function(name, expression) {
  file <- file.path(study_folder, paste0(name, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  value <- expression
  saveRDS(value, file)
  value
}


## =====================================================================
## Ground truths by adaptive integration (check 2.1)
## =====================================================================
## A different route per prior, none of them the package's:
##   flat_r          over g and over directions (the one-predictor prior as a
##                   mixture of normals, averaged over every direction)
##   flat_partial_r  over g: a g-prior mixture with g ~ Inverse-Gamma(1, n/2)
##   flat_r_squared  k >= 3: a g-prior mixture, g/n = B W, B ~ Beta(1, k/2 - 1),
##                   W ~ Inverse-Gamma(1, 1/2) (docs/appendix/A_derivations.md,
##                   section 8); k = 2: the same prior as flat_partial_r
##                   (a = k/2 = 1); k = 1: over the effect f directly
## Everything on the log scale, scaled by the largest value.

## The integral of exp(log_integrand(x)) over (lower, upper): the peak on a
## coarse grid, refined by optimize(); its width from the curvature; the
## interval split at the peak +- 1, 4, 10 and 30 widths, so that a narrow peak
## (very strong evidence) is never stepped over
log_integral_around_peak <- function(log_integrand, lower, upper, grid_step = 1) {
  grid   <- seq(lower, upper, by = grid_step)
  start  <- grid[which.max(log_integrand(grid))]
  peak   <- optimize(log_integrand, c(max(lower, start - grid_step), min(upper, start + grid_step)),
                     maximum = TRUE, tol = 1e-10)$maximum
  height <- log_integrand(peak)
  step   <- 1e-4 * max(1, abs(peak))
  curvature <- (log_integrand(peak + step) - 2 * height + log_integrand(peak - step)) / step^2
  width  <- if (is.finite(curvature) && curvature < 0) 1 / sqrt(-curvature) else grid_step
  breaks <- sort(unique(pmin(upper, pmax(lower, c(lower, upper,
                                                  peak + width * c(-30, -10, -4, -1, 0, 1, 4, 10, 30))))))
  scaled <- function(x) { value <- exp(log_integrand(x) - height); value[!is.finite(value)] <- 0; value }
  height + log(sum(vapply(seq_len(length(breaks) - 1), function(i)
    integrate(scaled, breaks[i], breaks[i + 1], rel.tol = 1e-10, subdivisions = 2000L,
              stop.on.error = FALSE)$value, numeric(1))))
}

## average over directions for one value of g: the angle phi between a
## direction and the data's direction, weight sin(phi)^(k - 2)
log_direction_average_by_integration <- function(scaled_share, sample_size, number_of_predictors) {
  half_n_minus_one <- (sample_size - 1) / 2
  log_value_at_zero <- -half_n_minus_one * log1p(-scaled_share)
  if (number_of_predictors == 1) return(log_value_at_zero)
  log_total_weight <- 0.5 * log(pi) - log(2) +
    lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
  scaled_integrand <- function(angle)
    exp(-half_n_minus_one * log1p(-scaled_share * cos(angle)^2) - log_value_at_zero +
        (number_of_predictors - 2) * log(sin(angle)) - log_total_weight)
  ## strong evidence makes a narrow peak near angle 0: split the interval there
  peak_width <- if (scaled_share > 0)
    sqrt((1 - scaled_share) / (2 * half_n_minus_one * scaled_share)) else Inf
  inner_breaks <- c(3, 10, 40) * peak_width
  break_points <- c(0, inner_breaks[inner_breaks < pi / 2], pi / 2)
  pieces <- vapply(seq_len(length(break_points) - 1), function(i)
    integrate(scaled_integrand, break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  log_value_at_zero + log(sum(pieces))
}

log_bayes_factor_flat_r_by_integration <- function(explained_share, sample_size, number_of_predictors) {
  if (number_of_predictors == 0) return(0)
  ## the integrand over log(g): IG(1, n/2) density times g, the g-prior factor,
  ## and the average over directions
  log_integrand <- function(log_g_values) vapply(log_g_values, function(log_g) {
    g <- exp(log_g)
    log(sample_size / 2) - 2 * log_g - sample_size / (2 * g) + log_g - 0.5 * log1p(g) +
      log_direction_average_by_integration(explained_share * g / (1 + g), sample_size,
                                           number_of_predictors)
  }, numeric(1))
  log_integral_around_peak(log_integrand, -40, 60)
}

## flat_partial_r: one integral over g, with g ~ Inverse-Gamma(1, n/2)
log_bf_g_inverse_gamma_by_integration <- function(r_squared, n, k) {
  log_integrand <- function(log_g) log(n / 2) - log_g - n / (2 * exp(log_g)) +
    (n - 1 - k) / 2 * log1p(exp(log_g)) - (n - 1) / 2 * log1p(exp(log_g) * (1 - r_squared))
  log_integral_around_peak(log_integrand, -40, 60)
}

## R-squared ~ Beta(a, 1) as a g-prior mixture, a < k/2: g/n = t = B W,
## B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2), p(t) = E[B exp(-B/(2t))] / (2 t^2)
log_hyperprior_density_g_over_n <- function(t, k, a) {
  weighted_beta <- function(b) b * dbeta(b, a, k / 2 - a) * exp(-b / (2 * t))
  split <- min(1, 60 * t)
  pieces <- integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
    (if (split < 1) integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
  log(pieces) - log(2) - 2 * log(t)
}
log_bf_g_mixture_by_integration <- function(r_squared, n, k, a) {
  log_integrand <- function(log_t) vapply(log_t, function(lt) {
    g <- n * exp(lt)
    value <- (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared)) +
      log_hyperprior_density_g_over_n(exp(lt), k, a) + lt
    if (is.finite(value)) value else -Inf
  }, numeric(1))
  log_integral_around_peak(log_integrand, -25, 25)
}

## log of the integral of v^m exp(-v^2/2 + b v) over v > 0, around its peak
log_I_by_integration <- function(b, m) {
  peak     <- (b + sqrt(b^2 + 4 * m)) / 2
  log_peak <- m * log(peak) - peak^2 / 2 + b * peak
  width    <- 1 / sqrt(1 + m / peak^2)           # of the integrand's peak
  scaled   <- function(v) exp(m * log(v) - v^2 / 2 + b * v - log_peak)
  breaks   <- c(max(0, peak - 40 * width), peak - width, peak + width, peak + 40 * width)
  breaks   <- sort(unique(pmax(breaks, 0)))
  log_peak + log(sum(vapply(seq_len(length(breaks) - 1), function(i)
    integrate(scaled, breaks[i], breaks[i + 1], rel.tol = 1e-10, stop.on.error = FALSE)$value,
    numeric(1))))
}

## flat_r_squared with one predictor, over the effect f: |f| has density
## 2 f / (1 + f^2)^2, and p(y | f) / p(y | 0) = exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0))
log_bf_one_predictor_flat_r_squared_by_integration <- function(r_squared, n) {
  m <- n - 2; r <- sqrt(r_squared)
  log_I0 <- log_I_by_integration(0, m)
  log_integrand <- function(f) vapply(f, function(f1) {
    b <- sqrt(n) * r * f1
    terms <- c(log_I_by_integration(b, m), log_I_by_integration(-b, m))
    -n * f1^2 / 2 + max(terms) + log(sum(exp(terms - max(terms)))) - log(2) - log_I0 +
      log(2 * f1) - 2 * log1p(f1^2)
  }, numeric(1))
  ## over log f, where the peak's width is of order 1 / sqrt(n) relative to f
  log_integral_around_peak(function(log_f) log_integrand(exp(log_f)) + log_f, log(1e-8), log(1e4),
                           grid_step = 0.25)
}

log_bf_by_integration <- function(prior, r_squared, n, k) {
  if (k == 0) return(0)
  switch(prior,
    flat_r         = log_bayes_factor_flat_r_by_integration(r_squared, n, k),
    flat_partial_r = log_bf_g_inverse_gamma_by_integration(r_squared, n, k),
    flat_r_squared = if (k == 1) log_bf_one_predictor_flat_r_squared_by_integration(r_squared, n) else
                     if (k == 2) log_bf_g_inverse_gamma_by_integration(r_squared, n, k) else
                     log_bf_g_mixture_by_integration(r_squared, n, k, a = 1))
}


## =====================================================================
## The likelihood ratio against its exact series, in 200-bit arithmetic
## =====================================================================
## All three priors average the same likelihood ratio (no prior in it) of an
## effect of length f in noise units against none, for q of k columns:
##   p(y | f) / p(y | 0) = exp(-n f^2 / 2) 1F1((m + 1)/2; q/2; n r^2 f^2 / 2),
##   m = n - k + q - 2, r^2 the columns' sample partial R^2
## (docs/appendix/A_derivations.md). The series of 1F1 has positive terms, so
## summing it in 200-bit arithmetic (Rmpfr) far past its largest term gives the
## exact value; the package computes the same ratio by an integral over the
## noise (R/likelihood.R in flatcor).
log_likelihood_ratio_by_series <- function(f, partial_r_squared, n, k, q, bits = 200) {
  a <- (n - k + q - 1) / 2; c_shape <- q / 2; z <- n * partial_r_squared * f^2 / 2
  largest_term <- (z + sqrt(z^2 + 4 * a * z)) / 2          # where the terms stop growing
  terms_needed <- ceiling(largest_term + 40 * sqrt(largest_term + 1) + 200)
  j <- Rmpfr::mpfr(0:terms_needed, bits)
  a_m <- Rmpfr::mpfr(a, bits); c_m <- Rmpfr::mpfr(c_shape, bits); z_m <- Rmpfr::mpfr(z, bits)
  log_terms <- lgamma(a_m + j) - lgamma(a_m) - lgamma(c_m + j) + lgamma(c_m) - lgamma(j + 1) +
    j * log(z_m)
  largest <- max(log_terms)
  -n * f^2 / 2 + as.numeric(largest + log(sum(exp(log_terms - largest))))
}


## =====================================================================
## Brute force, run until its Monte Carlo error is below 0.01 (check 2.3)
## =====================================================================
## The log Bayes factor without the formula: draw coefficient sets from the
## prior, see how well each predicts the outcome (the noise level averaged on
## a grid of 1 / noise), and average. Importance sampling from the flat_r
## recipe with weights 2a R^(2a - 1) for the other priors (bounded by 2a).
## Returns the per-draw log evidence, so the Monte Carlo error can be
## estimated: for log(mean(w)), MCSE = sd(w) / (mean(w) sqrt(N)).
log_evidence_per_draw_by_sampling <- function(predictors, outcome, prior, number_of_draws,
                                              draws_per_batch = 2000) {
  sample_size            <- length(outcome)
  centered_predictors    <- scale(predictors, scale = FALSE)
  centered_outcome       <- outcome - mean(outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  predictor_covariance   <- crossprod(centered_predictors) / sample_size
  inverse_noise_values <- seq(1e-4, 6 * sqrt(sample_size / outcome_sum_of_squares), length.out = 800)
  inverse_noise_step   <- inverse_noise_values[2] - inverse_noise_values[1]
  log_inverse_noise    <- log(inverse_noise_values)
  log_evidence_null <- log_sum_exp((sample_size - 2) * log_inverse_noise -
                                     inverse_noise_values^2 * outcome_sum_of_squares / 2) +
                       log(inverse_noise_step)
  a <- prior_shape(prior, ncol(predictors))
  unlist(lapply(seq_len(ceiling(number_of_draws / draws_per_batch)), function(batch) {
    coefficients    <- draw_coefficients_from_prior(draws_per_batch, predictor_covariance, prior = "flat_r")
    f_squared_drawn <- rowSums((coefficients %*% predictor_covariance) * coefficients)
    log_weight <- log(2 * a) + (2 * a - 1) / 2 * log(f_squared_drawn / (1 + f_squared_drawn))
    predicted_values <- centered_predictors %*% t(coefficients)
    overlap_with_outcome     <- drop(crossprod(centered_outcome, predicted_values))
    predicted_sum_of_squares <- colSums(predicted_values^2)
    log_integrand <- (sample_size - 2) *
        matrix(log_inverse_noise, draws_per_batch, length(inverse_noise_values), byrow = TRUE) -
      (outer(rep(outcome_sum_of_squares, draws_per_batch), inverse_noise_values^2) -
       2 * outer(overlap_with_outcome, inverse_noise_values) + predicted_sum_of_squares) / 2
    largest_per_row <- apply(log_integrand, 1, max)
    largest_per_row + log(rowSums(exp(log_integrand - largest_per_row))) + log(inverse_noise_step) +
      log_weight - log_evidence_null
  }))
}
estimate_from_log_evidence <- function(log_evidence) {
  largest <- max(log_evidence)
  weights <- exp(log_evidence - largest)
  c(log_bf = largest + log(mean(weights)), mcse = sd(weights) / (mean(weights) * sqrt(length(weights))))
}
## draws in rounds of 100,000 (in parallel) until the MCSE is below 0.01, at most
## 6.4 million; a smoke run takes one round of 20,000
log_bayes_factor_by_sampling_to_precision <- function(predictors, outcome, prior, target = 0.01) {
  round_size <- if (smoke_mode) 20000 else 100000
  maximum <- if (smoke_mode) round_size else 6.4e6
  log_evidence <- numeric(0)
  repeat {
    log_evidence <- c(log_evidence, unlist(parallel_map(seq_len(number_of_cores), function(core)
      log_evidence_per_draw_by_sampling(predictors, outcome, prior,
                                        ceiling(round_size / number_of_cores)))))
    estimate <- estimate_from_log_evidence(log_evidence)
    if (estimate[["mcse"]] < target || length(log_evidence) >= maximum) break
  }
  c(estimate, draws = length(log_evidence))
}


## =====================================================================
## Bridge sampling on the package's Stan fit (check 2.4)
## =====================================================================
## The joint-model Bayes factor against the intercept-only model as a ratio
## of marginal likelihoods. The full model's posterior draws come from the
## package's Stan model (system.file("stan", "flatcor.stan", package =
## "flatcor"), fitted with cmdstanr). bridgesampling's matrix method then
## needs the log posterior with every constant, which is written here in R
## (the Stan program drops constants):
##   likelihood    normal, with all constants
##   alpha, sigma  flat and 1/sigma (the package's improper priors; the same in
##                 both models, so their arbitrary scale cancels in the ratio)
##   R_magnitude   the size of R, density 2a R^(2a - 1) on (0, 1); with k = 1 the
##                 sign of R is + or - with probability 1/2, and the Stan model
##                 (flatcor 1.0.1) marginalises it: the likelihood is averaged
##                 over both signs, here as in the Stan program
##   z             k > 1: the auxiliary vector of the direction, density
##                 exp(-((|z| - 1)/0.1)^2 / 2) / C_k on R^k, with C_k = surface
##                 area of the unit sphere times the integral over the radius
##   coefficients  b = sigma f U^-1 (z/|z|), f = R / sqrt(1 - R^2), U the upper
##                 Cholesky factor of the predictors' sample covariance (as in
##                 the Stan model)
## The intercept-only model's marginal likelihood is in closed form:
##   n^(-1/2) (2 pi)^(-(n-1)/2) Gamma((n-1)/2) / 2 * (SST/2)^(-(n-1)/2).
## Its error: bridgesampling's coefficient of variation of the marginal
## likelihood, which is the MCSE of its log. Draws are doubled until it is
## below 0.01 (at most 4 x 32,000).
log_normalizing_constant_of_direction <- function(k) {
  radial <- integrate(function(r) exp(-((r - 1) / 0.1)^2 / 2 + (k - 1) * log(r)), 0, Inf,
                      rel.tol = 1e-12)$value
  log(2) + k / 2 * log(pi) - lgamma(k / 2) + log(radial)
}
log_marginal_likelihood_intercept_only <- function(outcome) {
  n <- length(outcome); sst <- sum((outcome - mean(outcome))^2)
  -0.5 * log(n) - (n - 1) / 2 * log(2 * pi) + lgamma((n - 1) / 2) - log(2) - (n - 1) / 2 * log(sst / 2)
}
bridge_sampling_log_bf <- function(predictors, outcome, prior, stan_model, seed, target = 0.01) {
  n <- length(outcome); k <- ncol(predictors); a <- prior_shape(prior, k)
  centered <- scale(predictors, scale = FALSE)
  upper_root <- chol(crossprod(centered) / n)
  log_c <- if (k > 1) log_normalizing_constant_of_direction(k) else 0
  parameter_names <- c("alpha", "sigma", "R_magnitude", if (k > 1) paste0("z[", seq_len(k), "]"))
  log_posterior <- function(pars, data) {
    sigma <- pars[["sigma"]]; R <- pars[["R_magnitude"]]
    direction <- if (k > 1) { z <- pars[paste0("z[", seq_len(k), "]")]; z / sqrt(sum(z^2)) } else 1
    b <- sigma * R / sqrt(1 - R^2) * backsolve(upper_root, direction)
    log_prior <- -log(sigma) + log(2 * a) + (2 * a - 1) * log(R) +
      (if (k > 1) -((sqrt(sum(z^2)) - 1) / 0.1)^2 / 2 - log_c else 0)
    log_likelihood <- function(sign) sum(dnorm(outcome, pars[["alpha"]] + sign * drop(centered %*% b), sigma, log = TRUE))
    if (k > 1) return(log_likelihood(1) + log_prior)
    log_sum_exp(c(log_likelihood(1), log_likelihood(-1))) - log(2) + log_prior
  }
  lower <- setNames(c(-Inf, 0, 0, rep(-Inf, if (k > 1) k else 0)), parameter_names)
  upper <- setNames(c(Inf, Inf, 1, rep(Inf, if (k > 1) k else 0)), parameter_names)
  draws_per_chain <- if (smoke_mode) 1000 else 4000
  repeat {
    fit <- stan_model$sample(data = list(N = n, K = k, X = predictors, y = outcome, prior_only = 0L,
                                         proper_stand_ins = 0L, r_squared_shape = a),
                             chains = 4, parallel_chains = min(4, number_of_cores),
                             iter_warmup = 1000, iter_sampling = draws_per_chain,
                             adapt_delta = 0.99, seed = seed, refresh = 0, show_messages = FALSE)
    samples <- posterior::as_draws_matrix(fit$draws(parameter_names))
    samples <- unclass(samples)[, parameter_names, drop = FALSE]
    bridge <- suppressWarnings(bridgesampling::bridge_sampler(samples = samples,
                log_posterior = log_posterior, data = NULL, lb = lower, ub = upper, silent = TRUE))
    error <- bridgesampling::error_measures(bridge)$cv
    if (error < target || smoke_mode || draws_per_chain >= 32000) break
    draws_per_chain <- 2 * draws_per_chain
  }
  c(log_bf = bridge$logml - log_marginal_likelihood_intercept_only(outcome), error = error,
    draws = 4 * draws_per_chain, divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
}


## =====================================================================
## Data: the example and two more datasets, for the routes that need data
## =====================================================================
set.seed(42)
example_n <- 60
example_correlations <- matrix(0.4, 3, 3); diag(example_correlations) <- 1
example_predictors <- matrix(rnorm(example_n * 3), example_n) %*% chol(example_correlations)
colnames(example_predictors) <- c("x1", "x2", "x3")
example_outcome <- drop(example_predictors %*% c(0.3, 0.2, 0)) + rnorm(example_n)
example_data <- data.frame(example_predictors, y = example_outcome)

data_cases <- list(
  list(label = "k = 3, correlated (example data)", predictors = example_predictors,
       outcome = example_outcome),
  local({ set.seed(11); x <- matrix(rnorm(80), 80, 1)
          list(label = "k = 1, n = 80", predictors = x, outcome = drop(0.3 * x) + rnorm(80)) }),
  local({ set.seed(12); x <- matrix(rnorm(100 * 5), 100) %*% chol(0.3 + 0.7 * diag(5))
          list(label = "k = 5, n = 100, correlated", predictors = x,
               outcome = drop(x %*% c(0.2, 0, 0.1, 0, 0)) + rnorm(100)) }))
if (smoke_mode) data_cases <- data_cases[1]


## =====================================================================
## THE CHECKS
## =====================================================================
sweep_settings <- subset(expand.grid(explained_share      = c(0.001, 0.01, 0.1, 0.3, 0.6, 0.9, 0.99, 0.9999),
                                     sample_size          = c(10, 30, 100, 1000, 10000, 1e5),
                                     number_of_predictors = c(1, 2, 3, 5, 10, 30)),
                         sample_size > number_of_predictors + 1)
if (smoke_mode) sweep_settings <- subset(sweep_settings, sample_size %in% c(30, 1000) &
                                           number_of_predictors %in% c(1, 3) &
                                           explained_share %in% c(0.01, 0.3, 0.9))

stan_model <- local({
  ## compiled from a copy of the package's Stan file in the cache, so the
  ## executable is not written into the installed package
  stan_file <- file.path(study_folder, "flatcor.stan")
  file.copy(system.file("stan", "flatcor.stan", package = "flatcor"), stan_file, overwrite = TRUE)
  cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
})

sweeps <- list(); brute_results <- list(); bridge_results <- list()
for (prior in priors) {
  cat("\n## ========================= ", prior, " =========================\n\n", sep = "")

  ## 2.1 adaptive integration of the g-prior form
  sweep <- cached(paste0("sweep_", prior), {
    message(prior, ": integrate() for ", nrow(sweep_settings), " settings ...")
    values <- parallel_map(seq_len(nrow(sweep_settings)), function(row) {
      s <- sweep_settings[row, ]
      c(code         = log_bf(prior, s$explained_share, s$sample_size, s$number_of_predictors),
        ground_truth = log_bf_by_integration(prior, s$explained_share, s$sample_size,
                                             s$number_of_predictors))
    })
    cbind(sweep_settings, do.call(rbind, values))
  })
  sweep$error <- abs(sweep$code - sweep$ground_truth)
  sweep$relative_error <- sweep$error / pmax(1, abs(sweep$ground_truth))
  sweeps[[prior]] <- transform(sweep, prior = prior)
  ordinary <- abs(sweep$ground_truth) < 50
  record_check(sprintf("2.1a %s: = integrate() of the g-prior form, ordinary evidence (%d settings)",
                       prior, sum(ordinary)),
               max(sweep$error[ordinary]) < 1e-5, sprintf("largest error %.1e", max(sweep$error[ordinary])))
  worst <- which.max(sweep$relative_error)
  record_check(sprintf("2.1b %s: = integrate() of the g-prior form, all settings, relative (%d)",
                       prior, nrow(sweep)),
               max(sweep$relative_error) < 1e-5 &&
                 all(sweep$relative_error[abs(sweep$ground_truth) > 1] < 1e-6 |
                       sweep$error[abs(sweep$ground_truth) > 1] < 1e-5),
               sprintf("largest relative error %.1e (R^2 %g, n %g, k %g, log BF %.3g)",
                       sweep$relative_error[worst], sweep$explained_share[worst], sweep$sample_size[worst],
                       sweep$number_of_predictors[worst], sweep$ground_truth[worst]))

  ## 2.3 brute force, to MCSE 0.01
  brute <- cached(paste0("brute_force_", prior), {
    RNGkind("L'Ecuyer-CMRG"); set.seed(21)
    lapply(data_cases, function(case) {
      message(prior, ": brute force, ", case$label, " ...")
      log_bayes_factor_by_sampling_to_precision(case$predictors, case$outcome, prior)
    })
  })
  for (i in seq_along(data_cases)) {
    case <- data_cases[[i]]; estimate <- brute[[i]]
    package_value <- log_bf(prior, sample_r_squared(case$predictors, case$outcome),
                            length(case$outcome), ncol(case$predictors))
    brute_results[[length(brute_results) + 1]] <- data.frame(prior = prior, data = case$label,
      package = package_value, brute_force = estimate[["log_bf"]], mcse = estimate[["mcse"]])
    record_check(sprintf("2.3 %s: brute force, %s", prior, case$label),
                 abs(package_value - estimate[["log_bf"]]) < 3 * estimate[["mcse"]] &&
                   (smoke_mode || estimate[["mcse"]] < 0.01),
                 sprintf("package %.4f, brute force %.4f (MCSE %.4f, %s draws)", package_value,
                         estimate[["log_bf"]], estimate[["mcse"]], format(estimate[["draws"]], big.mark = ",")))
  }

  ## 2.4 bridge sampling on the package's Stan fit, to an error of 0.01
  bridged <- cached(paste0("bridge_sampling_", prior), lapply(seq_along(data_cases), function(i) {
    message(prior, ": bridge sampling, ", data_cases[[i]]$label, " ...")
    bridge_sampling_log_bf(data_cases[[i]]$predictors, data_cases[[i]]$outcome, prior, stan_model,
                           seed = 600 + i)
  }))
  for (i in seq_along(data_cases)) {
    case <- data_cases[[i]]; estimate <- bridged[[i]]
    package_value <- log_bf(prior, sample_r_squared(case$predictors, case$outcome),
                            length(case$outcome), ncol(case$predictors))
    bridge_results[[length(bridge_results) + 1]] <- data.frame(prior = prior, data = case$label,
      package = package_value, bridge_sampling = estimate[["log_bf"]], error = estimate[["error"]],
      divergent = estimate[["divergent"]])
    record_check(sprintf("2.4 %s: bridge sampling on the Stan fit, %s", prior, case$label),
                 abs(package_value - estimate[["log_bf"]]) < 3 * estimate[["error"]] &&
                   (smoke_mode || estimate[["error"]] < 0.01),
                 sprintf("package %.4f, bridge sampling %.4f (error %.4f, %s draws, %d divergent)",
                         package_value, estimate[["log_bf"]], estimate[["error"]],
                         format(estimate[["draws"]], big.mark = ","), as.integer(estimate[["divergent"]])))
  }
}

## 2.2 the likelihood ratio the three priors share, against its 200-bit series
cat("\n## ========================= the shared likelihood ratio =========================\n\n")
series_settings <- subset(expand.grid(f = c(0.05, 0.2, 0.5, 1), partial_r_squared = c(0.01, 0.3, 0.9),
                                      sample_size = c(10, 100, 1000, 10000), number_of_predictors = c(1, 3, 10, 30),
                                      columns = c("one", "all"), stringsAsFactors = FALSE),
                          sample_size > number_of_predictors + 2 &
                            sample_size * partial_r_squared * f^2 / 2 < 2e4)
if (smoke_mode) series_settings <- subset(series_settings, sample_size == 100 & number_of_predictors %in% c(1, 3))
series <- cached("likelihood_ratio_series", {
  message("200-bit series for ", nrow(series_settings), " settings ...")
  values <- parallel_map(seq_len(nrow(series_settings)), function(row) {
    s <- series_settings[row, ]
    q <- if (s$columns == "one") 1 else s$number_of_predictors
    c(code   = log_likelihood_ratio_block(s$f, s$partial_r_squared, s$sample_size, s$number_of_predictors, q),
      series = log_likelihood_ratio_by_series(s$f, s$partial_r_squared, s$sample_size,
                                              s$number_of_predictors, q))
  })
  cbind(series_settings, do.call(rbind, values))
})
series$error <- abs(series$code - series$series)
record_check(sprintf("2.2 the shared likelihood ratio = its 200-bit series (%d settings, q = 1 and q = k)",
                     nrow(series)),
             max(series$error) < 1e-6, sprintf("largest error %.1e", max(series$error)))


save_study_results("02", list(sweep = do.call(rbind, sweeps), likelihood_ratio = series,
                              brute_force = do.call(rbind, brute_results),
                              bridge_sampling = do.call(rbind, bridge_results)))
report_checks()
