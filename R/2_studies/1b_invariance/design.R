## =====================================================================
##  Study 1b: is the Bayes factor invariant to units and coding?   (heading 1, row b)
##  The design. What is tested, what is expected and why:
##  docs/validation_tests.md, study 1b. The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (3) × dataset (ordinary regression, R/1_data/2_ordinary_regression.R) ×
##   change (rescaled, recoded)
priors_1b <- data.frame(prior = c("flat_partial_r", "flat_r", "flat_r_squared"))


## ---- How many (docs/validation_tests.md, "How many datasets") -------------------------
## Δ per quantity: the Monte Carlo SE must be at most Δ / 1.96
deltas_1b <- c(
  "change in log BF, rescaled" = log(1.1),   # .095: the Bayes factor off by 10%
  "change in log BF, recoded"  = log(1.1)
)
## The first round (the Chow-Robbins rule's, per prior); later rounds add
## datasets until both changes reach their targets (R/shared/rounds.R)
first_round_1b <- datasets(50, in_smoke = 2)


## ---- Bookkeeping (no effect on the results) ------------------------------------
batch_size_1b <- if (smoke) 2 else 5   # datasets per piece of work, saved as it finishes
