## =====================================================================
##  Study 1a: does the software draw the prior it claims?   (heading 1, row a)
##  The design. What is tested, what is expected and why, and the reasons
##  for each choice: docs/validation_tests.md, study 1a. The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (flat R, flat partial r, flat R²) × predictor design (one, three
##   correlated, eight; R/1_data/1_prior_only_designs.R) = 9 settings
settings_1a <- expand.grid(
  prior = c("flat_r", "flat_partial_r", "flat_r_squared"),
  design = c("one", "three_correlated", "eight"),
  stringsAsFactors = FALSE
)
percentiles_1a <- c(0.05, 0.25, 0.5, 0.75, 0.95)


## ---- How many (docs/validation_tests.md, "How many datasets") -------------------------
## Δ for every percentile: 10% of the quantity's SD under the prior (its
## natural unit; also for a median of 0, as the partial correlations have).
## The chains are lengthened until every percentile's Monte Carlo SE is at
## most Δ / 1.96, starting from start_draws_1a draws in all.
delta_share_1a    <- 0.1
start_draws_1a    <- if (smoke) 500 else 10000
reference_size_1a <- if (smoke) 20000 else 200000   # draws of the prior as a recipe (2nd route)
