## =====================================================================
##  Study 1b: is the Bayes factor invariant to units and coding?   (heading 1, row b)
##  The design. What is tested, what is expected and why:
##  docs/validation_tests.md, study 1b. The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (3) × dataset (ordinary regression, R/1_data/2_ordinary_regression.R) ×
##   change (rescaled, recoded)
priors_1b <- data.frame(prior = c("flat_partial_r", "flat_r", "flat_r_squared"))


## ---- How many --------------------------------------------------------------------
datasets_1b <- datasets(20, in_smoke = 2)   # per prior (set by the target)


## ---- Bookkeeping (no effect on the results) ------------------------------------
batch_size_1b <- if (smoke) 2 else 5   # datasets per piece of work, saved as it finishes
