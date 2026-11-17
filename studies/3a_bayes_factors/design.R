## =====================================================================
##  Study 3a: are the Bayes factors right?   (heading 3, row a)
##  The design. What is tested, what is expected and why, how a dataset is
##  drawn and the reasons for each choice: docs/validation_tests.md, study 3a.
##  The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (3) × setting (k = 3, n = 50; k = 8, n = 200) = 6 cells,
##   data from the prior (R/data/from_the_prior.R), half from the model
##   without x1
cells_3a <- data.frame(
  prior = rep(c("flat_partial_r", "flat_r", "flat_r_squared"), each = 2),
  k     = c(3, 8),
  n     = c(50, 200)
)
if (smoke) cells_3a <- cells_3a[cells_3a$k == 3, ]


## ---- How many --------------------------------------------------------------------
datasets_per_cell_3a <- datasets(1000, in_smoke = 40)   # part A (set by cost; provisional)
repeat_datasets_3a   <- datasets(20, in_smoke = 2)      # part B (set by the target)
repeat_runs_3a       <- 5                               # part B: runs per dataset


## ---- Bookkeeping (no effect on the results) ------------------------------------
batch_size_3a <- if (smoke) 4 else 25   # datasets per piece of work, saved as it finishes
