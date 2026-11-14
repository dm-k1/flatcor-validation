## =====================================================================
##  Study 1: is the prior sampled correctly?
##  (docs/appendix/C_study_design.md, study 1)
##
##  Aims. flatcor's Stan model (inst/stan/flatcor.stan in the installed
##  package, used by fit_regression()) samples R, with R-squared ~ Beta(a, 1),
##  and a uniform direction, and computes the coefficients from them. Does it
##  reproduce each of the three priors?
##  Data. None: the model samples the prior only (prior_only = 1), with proper
##  stand-in priors on the intercept and the noise. The predictors are fixed
##  designs (n = 60): one predictor, three correlated predictors
##  (correlation .4) and eight uncorrelated ones, each drawn once from its own
##  stream (R/data/data_sets.R).
##  Estimands. The distribution of R-squared, of the correlation (k = 1) and of
##  each coefficient's partial correlation, against its exact form.
##  Methods. The Stan model for flat_r (a = 1/2), flat_partial_r (a = k/2) and
##  flat_r_squared (a = 1).
##  Performance. The Kolmogorov-Smirnov distance between the draws and the
##  exact distribution, judged against the KS critical value at the .01 level
##  for the effective number of draws, 1.63 / sqrt(ESS) (two-sample: for the
##  ESS and the size of the reference sample). The ESS is posterior::ess_bulk()
##  of the checked quantity, so autocorrelated draws count for what they are
##  worth.
##  n_sim. Enough iterations that every checked quantity reaches an ESS of
##  n_sim_of("ks_draws") = 26,569 (critical distance <= .01;
##  R/shared/precision.R): the chains are lengthened until it does.
##
##  Checks, the same four for each prior (a its shape: flat_r 1/2,
##  flat_partial_r k/2, flat_r_squared 1):
##    1.1 k = 1: r^2 ~ Beta(a, 1), a for k = 1 (r^2 rather than r: under
##        flat_r_squared r is bimodal near -1 and 1, so its sign mixes slowly)
##    1.2 k = 3 correlated: R-squared ~ Beta(a, 1)
##    1.3 k = 8: R-squared ~ Beta(a, 1)
##    1.4 k = 3 correlated: each partial correlation matches the prior as a
##        recipe, draw_coefficients_from_prior() (two-sample)
##  Each check passes when the KS distance is below the critical value for
##  its ESS.

##  Run from the project folder:  Rscript R/studies/01_prior_sampling.R
##  Smoke test (ESS target 500, about a minute):
##    VALIDATION_SMOKE=1 Rscript R/studies/01_prior_sampling.R
##  Raw results are cached in cache/study_01/flatcor_<build>/ (delete it to rerun); the
##  tables and checks are saved to results/studies/study_01.rds.
## =====================================================================
library(posterior)
source("R/shared/settings.R")             # flatcor, cores, record_check()
source("R/shared/precision.R")            # n_sim_of(), smoke_mode
source("R/data/data_sets.R")            # with_stream_of()

study_folder <- study_cache_folder("01")

## the ESS every checked quantity must reach (smoke mode: 500)
target_ess <- if (smoke_mode) 500L else precision_targets$ks_draws$n_sim
chains <- 4


## ---- Sampling helper: cmdstanr if installed, otherwise rstan -------------------
## Returns the draws of b and sigma as a posterior::draws_array (iterations x
## chains x variables), so the ESS of any derived quantity can be computed.
stan_model_code <- paste(readLines(system.file("stan", "flatcor.stan", package = "flatcor")),
                         collapse = "\n")
use_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
prior_model <- if (use_cmdstanr) {
  cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_model_code))
} else {
  rstan::stan_model(model_code = stan_model_code)
}
sample_prior_draws <- function(stan_data, iter_sampling, seed) {
  if (use_cmdstanr) {
    fit <- prior_model$sample(data = stan_data, chains = chains, parallel_chains = min(chains, number_of_cores),
                              iter_warmup = 1000, iter_sampling = iter_sampling, seed = seed,
                              refresh = 0, show_messages = FALSE)
    fit$draws(c("b", "sigma"))
  } else {
    fit <- rstan::sampling(prior_model, data = stan_data, chains = chains,
                           iter = 1000 + iter_sampling, warmup = 1000, seed = seed, refresh = 0)
    as_draws_array(as.array(fit, pars = c("b", "sigma")))
  }
}


## ---- The fixed designs ----------------------------------------------------------
sample_size <- 60
covariance_of <- function(predictors) {
  centered <- scale(predictors, scale = FALSE)
  crossprod(centered) / nrow(centered)
}
designs <- list(
  one = with_stream_of("study01-design-one", matrix(rnorm(sample_size), sample_size)),
  three_correlated = with_stream_of("study01-design-three-correlated",
    matrix(rnorm(sample_size * 3), sample_size) %*%
      chol(matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), 3, 3))),
  eight = with_stream_of("study01-design-eight", matrix(rnorm(sample_size * 8), sample_size)))


## ---- Quantities computed from the draws (each an iterations x chains matrix) ----
coefficient_draws <- function(draws, j) as.matrix(extract_variable_array(draws, sprintf("b[%d]", j))[, , 1])
sigma_draws <- function(draws) as.matrix(extract_variable_array(draws, "sigma")[, , 1])

## explained share: f^2 = b' Sigma b / sigma^2, R^2 = f^2 / (1 + f^2)
implied_r_squared <- function(draws, predictor_covariance) {
  k <- ncol(predictor_covariance)
  coefficients <- sapply(seq_len(k), function(j) as.vector(coefficient_draws(draws, j)))
  f_squared <- rowSums((coefficients %*% predictor_covariance) * coefficients) /
               as.vector(sigma_draws(draws))^2
  matrix(f_squared / (1 + f_squared), ncol = chains)
}

## each coefficient's partial correlation: the coefficient in noise units times
## the predictor's spread left after the others, rho = f_j / sqrt(1 + f_j^2)
partial_correlation <- function(draws, predictor_covariance, j) {
  spread_left <- sqrt(1 / diag(solve(predictor_covariance)))[j]
  f_j <- coefficient_draws(draws, j) / sigma_draws(draws) * spread_left
  f_j / sqrt(1 + f_j^2)
}

## the signed correlation with one predictor
correlation_one <- function(draws, predictor_covariance) {
  f <- coefficient_draws(draws, 1) * sqrt(predictor_covariance[1, 1]) / sigma_draws(draws)
  f / sqrt(1 + f^2)
}


## ---- Sampling until every checked quantity reaches the target ESS --------------
## Each setting (design and prior) is sampled with a fixed seed. If the
## smallest ESS of its checked quantities falls short of the target, the
## chains are lengthened in proportion and sampled again (same seed); the rule
## concerns precision only, not the outcome of any check.
priors <- c("flat_r", "flat_partial_r", "flat_r_squared")
settings <- unlist(lapply(seq_along(priors), function(i) lapply(seq_along(designs), function(j)
  list(name = paste(priors[i], names(designs)[j], sep = "-"), design = names(designs)[j],
       prior = priors[i], seed = 100 + 10 * i + j))), recursive = FALSE)

## the quantities whose ESS must reach the target, per setting
checked_quantities <- function(draws, setting) {
  covariance <- covariance_of(designs[[setting$design]])
  k <- ncol(covariance)
  ## k = 1: r^2 (under flat_r_squared r is bimodal near -1 and 1, so its
  ## sign mixes slowly while r^2 does not)
  if (k == 1) return(list(correlation_one(draws, covariance)^2))
  c(list(implied_r_squared(draws, covariance)),
    if (setting$design == "three_correlated")
      lapply(seq_len(k), function(j) partial_correlation(draws, covariance, j)))
}

sample_setting <- function(setting) {
  cache_file <- file.path(study_folder, paste0(setting$name, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  predictors <- designs[[setting$design]]
  stan_data <- list(N = sample_size, K = ncol(predictors), X = predictors,
                    y = with_stream_of(paste0("study01-outcome-", setting$design), rnorm(sample_size)),
                    prior_only = 1L, proper_stand_ins = 1L,
                    r_squared_shape = prior_shape(setting$prior, ncol(predictors)))
  iter_sampling <- ceiling(1.2 * target_ess / chains)
  for (attempt in 1:4) {
    draws <- sample_prior_draws(stan_data, iter_sampling, setting$seed)
    smallest_ess <- min(vapply(checked_quantities(draws, setting), ess_bulk, numeric(1)))
    if (smallest_ess >= target_ess) break
    iter_sampling <- ceiling(iter_sampling * 1.2 * target_ess / smallest_ess)
  }
  result <- list(draws = draws, iter_sampling = iter_sampling, smallest_ess = smallest_ess)
  saveRDS(result, cache_file)
  result
}

samples <- setNames(lapply(settings, sample_setting), vapply(settings, `[[`, "", "name"))
cat("\nIterations per chain and smallest ESS per setting (target", target_ess, "):\n")
print(data.frame(setting = names(samples),
                 iterations_per_chain = vapply(samples, `[[`, numeric(1), "iter_sampling"),
                              smallest_ess = round(vapply(samples, `[[`, numeric(1), "smallest_ess"))),
      row.names = FALSE)


## ---- The KS check against the critical value for the ESS ------------------------
## one-sample: critical distance 1.63 / sqrt(ESS); two-sample against a
## reference sample of size m: 1.63 * sqrt((ESS + m) / (ESS * m))
ks_against <- function(quantity, reference, reference_size = NULL) {
  distance <- suppressWarnings(unname(ks.test(as.vector(quantity), reference)$statistic))
  ess <- ess_bulk(quantity)
  critical <- if (is.null(reference_size)) 1.63 / sqrt(ess) else
    1.63 * sqrt((ess + reference_size) / (ess * reference_size))
  c(distance = distance, ess = ess, critical = critical)
}
record_ks <- function(label, results) {
  results <- if (is.list(results)) results else list(results)
  distances <- vapply(results, `[[`, numeric(1), "distance")
  criticals <- vapply(results, `[[`, numeric(1), "critical")
  record_check(label, all(distances < criticals),
               sprintf("KS %s; critical %s (ESS %s)",
                       paste(sprintf("%.4f", distances), collapse = " "),
                       paste(sprintf("%.4f", criticals), collapse = " "),
                       paste(round(vapply(results, `[[`, numeric(1), "ess")), collapse = " ")))
}
covariance <- lapply(designs, covariance_of)

reference_size <- if (smoke_mode) 20000 else 200000
spread_left <- sqrt(1 / diag(solve(covariance$three_correlated)))
for (prior in priors) {
  shape_of <- function(k) prior_shape(prior, k)
  draws_of <- function(design) samples[[paste(prior, design, sep = "-")]]$draws

  ## 1.1 k = 1: r^2 ~ Beta(a, 1)
  record_ks(sprintf("1.1 %s, k = 1: r^2 ~ Beta(%g, 1)", prior, shape_of(1)),
            ks_against(correlation_one(draws_of("one"), covariance$one)^2,
                       function(r2) pbeta(r2, shape_of(1), 1)))

  ## 1.2, 1.3 R-squared ~ Beta(a, 1) at k = 3 (correlated) and k = 8
  record_ks(sprintf("1.2 %s, k = 3 correlated: R-squared ~ Beta(%g, 1)", prior, shape_of(3)),
            ks_against(implied_r_squared(draws_of("three_correlated"), covariance$three_correlated),
                       function(r2) pbeta(r2, shape_of(3), 1)))
  record_ks(sprintf("1.3 %s, k = 8: R-squared ~ Beta(%g, 1)", prior, shape_of(8)),
            ks_against(implied_r_squared(draws_of("eight"), covariance$eight),
                       function(r2) pbeta(r2, shape_of(8), 1)))

  ## 1.4 k = 3 correlated: each partial correlation matches the prior as a
  ## recipe (draw_coefficients_from_prior(), in noise units)
  recipe <- with_stream_of(paste0("study01-recipe-", prior, "-three"),
    draw_coefficients_from_prior(reference_size, covariance$three_correlated, prior = prior))
  recipe_partials <- sweep(recipe, 2, spread_left, "*")
  recipe_partials <- recipe_partials / sqrt(1 + recipe_partials^2)
  record_ks(sprintf("1.4 %s, k = 3 correlated: each partial correlation matches the recipe", prior),
            lapply(1:3, function(j)
              ks_against(partial_correlation(draws_of("three_correlated"), covariance$three_correlated, j),
                         recipe_partials[, j], reference_size = reference_size)))
}

save_study_results("01", list(iterations = vapply(samples, `[[`, numeric(1), "iter_sampling"),
                              smallest_ess = vapply(samples, `[[`, numeric(1), "smallest_ess")))
report_checks()
