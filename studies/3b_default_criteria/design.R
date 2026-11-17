## =====================================================================
##  Study 3b: do the Bayes factors meet the default-prior criteria?   (heading 3, row b)
##  The design. What is tested (3b.1 ... 3b.9), what is expected and why:
##  docs/validation_tests.md, study 3b. The code: code.R.
## =====================================================================


priors_3b <- c("flat_partial_r", "flat_r", "flat_r_squared")


## ---- 3b.1 predictive matching, 3b.2 information consistency ---------------------
##   3b.1: prior (3) × k (1, 2, 3) at n = k + 1 = 9 cells
##   3b.2: prior (3) × R² (.9, .99, .999, .9999) at n = 20, k = 3 = 12 cells
##   one dataset each, with a set sample R-squared (R/data/set_r_squared.R)
formula_checks_3b <- rbind(
  expand.grid(
    item = "3b.1", k = 1:3, r_squared = 1, prior = priors_3b,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    item = "3b.2", k = 3, r_squared = c(0.9, 0.99, 0.999, 0.9999), prior = priors_3b,
    stringsAsFactors = FALSE
  )
)
formula_checks_3b$n <- ifelse(formula_checks_3b$item == "3b.1", formula_checks_3b$k + 1, 20)


## ---- 3b.3-3b.6 the rates ----------------------------------------------------------
##   setting (3b.3: one null next to two real; 3b.4: x1 real, q = 1, 3, 6;
##   3b.5: intercept only, q = 1, 3, 6; 3b.6: alternative) × n (50, 200,
##   1000, 5000; the alternative from 200) × prior (3) = 93 cells
##   (R/data/nested_models.R)
rate_settings_3b <- rbind(
  data.frame(item = "3b.3", setting = "x3 null next to two real", k_small = 2, q = 1),
  data.frame(item = "3b.4", setting = "x1 real",                  k_small = 1, q = c(1, 3, 6)),
  data.frame(item = "3b.5", setting = "intercept only",           k_small = 0, q = c(1, 3, 6)),
  data.frame(item = "3b.6", setting = "alternative",              k_small = 2, q = 1)
)

rate_cells_3b <- merge(rate_settings_3b, data.frame(n = c(50, 200, 1000, 5000)))
rate_cells_3b <- merge(rate_cells_3b, data.frame(prior = priors_3b, stringsAsFactors = FALSE))
rate_cells_3b <- subset(rate_cells_3b, setting != "alternative" | n >= 200)
if (smoke) {
  rate_cells_3b <- subset(rate_cells_3b, q == 1 & n %in% c(50, 200) & setting != "alternative")
}

rate_datasets_3b <- datasets(100, in_smoke = 2)   # per cell (set by cost; provisional)


## ---- 3b.7 intrinsic consistency (the prior as a recipe; no Bayes factors) --------
##   prior (3) × n (30, 300, 100,000) = 9 cells
intrinsic_n_3b     <- c(30, 300, 100000)
intrinsic_draws_3b <- if (smoke) 5000 else 26569


## ---- 3b.8-3b.9 k growing with n (k at most 90) --------------------------------------
##   growth (k = n / 4 at n = 50, 100, 200, 360; k = sqrt(n) at n = 50 ... 2,000)
##   × truth (null, a small model) × method (3 priors, JZS, Zellner-Siow)
##   = 100 cells (R/data/growing_k.R)
growing_methods_3b <- c("flat partial r", "flat R", "flat R-squared", "JZS default", "Zellner-Siow")

growing_cells_3b <- merge(growing_k_cells(), data.frame(truth = c("null", "true")))
growing_cells_3b <- merge(growing_cells_3b, data.frame(method = growing_methods_3b))
if (smoke) {
  growing_cells_3b <- subset(growing_cells_3b, n %in% c(50, 100) & growth == "k = n / 4")
}

growing_datasets_3b <- datasets(50, in_smoke = 2)   # per cell (set by cost; provisional)


## ---- For every part ----------------------------------------------------------------
bootstrap_replicates_3b <- if (smoke) 100 else 1000   # for the MCSE of medians and slopes


## ---- Bookkeeping (no effect on the results) ------------------------------------
batch_size_3b <- if (smoke) 2 else 10   # datasets per piece of work, saved as it finishes
