## =====================================================================
##  Study 2a: is the posterior computed correctly?   (heading 2, row a)
##  The design. What is tested, what is expected and why:
##  docs/validation_tests.md, study 2a. The code: code.R.
## =====================================================================


## ---- The design -------------------------------------------------------------------
##   prior (3) × setting (k = 8, n = 60; k = 3, n = 50) = 6 cells,
##   data from the prior (R/data/from_the_prior.R)
settings_2a <- merge(
  data.frame(prior = c("flat_partial_r", "flat_r", "flat_r_squared")),
  data.frame(k = c(8, 3), n = c(60, 50))
)
if (smoke) settings_2a <- settings_2a[settings_2a$k == 8, ]

iterations_2a <- if (smoke) {
  c(warmup = 300, sampling = 300)
} else {
  c(warmup = 1000, sampling = 1000)
}


## ---- How many --------------------------------------------------------------------
fits_2a <- datasets(1000, in_smoke = 10)   # fits per cell (set by the target)


## ---- Bookkeeping (no effect on the results) ------------------------------------
batch_size_2a <- if (smoke) 10 else 50   # fits per piece of work, saved as it finishes
