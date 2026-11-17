## =====================================================================
##  Study 1a: does the software draw the prior it claims?   (heading 1, row a)
##  The design. What is tested, what is expected and why, and the reasons
##  for each choice: docs/validation_tests.md, study 1a. The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (flat R, flat partial r, flat R²) × predictor design (one, three
##   correlated, eight; R/data/prior_only_designs.R) = 9 settings
settings_1a <- expand.grid(
  prior = c("flat_r", "flat_partial_r", "flat_r_squared"),
  design = c("one", "three_correlated", "eight"),
  stringsAsFactors = FALSE
)
percentiles_1a <- c(0.05, 0.25, 0.5, 0.75, 0.95)


## ---- How many --------------------------------------------------------------------
## effective draws per quantity (set by the target), and draws of the prior
## as a recipe (the second route)
target_ess_1a     <- if (smoke) 500 else 26569
reference_size_1a <- if (smoke) 20000 else 200000
