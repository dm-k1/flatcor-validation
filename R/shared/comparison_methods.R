## =====================================================================
##  The methods compared in the validation studies: one set, used
##  everywhere.
##
##  The three priors, all R-squared ~ Beta(a, 1) with a uniform direction
##  (flatcor):
##    flat partial r   a = k/2: every correlation matrix of (y, x1 ... xk)
##                     equally likely; rstanarm stan_lm's prior with eta = 1
##    flat R           a = 1/2: R flat for every k
##    flat R-squared   a = 1:   R-squared flat for every k
##  Two references, the standard defaults (BayesFactor / JASP; Liang et al., 2008):
##    JZS default      Cauchy prior, scale sqrt(2)/4 ("medium")
##    Zellner-Siow     Cauchy prior, scale 1 (R-squared ~ Beta(k/2, 1/2))
##
##  Every Bayes factor is a comparison of two models, each with its own prior:
##  for one coefficient (or one term), the model with it against the model
##  without it.
## =====================================================================

## ---- Method names and colours (the same in every study) --------------------
contenders <- c(`flat partial r` = "flat_partial_r", `flat R` = "flat_r",
                `flat R-squared` = "flat_r_squared")
model_methods <- c(names(contenders), "JZS default", "Zellner-Siow")

## Colour for the figures: one colour per prior, the references in grey.
## Names may carry a suffix ("JZS (lmBF)").
method_colour_of <- function(method)
  if (startsWith(method, "flat partial r")) "#0072B2" else
  if (startsWith(method, "flat R-squared")) "#009E73" else
  if (startsWith(method, "flat R"))         "#D55E00" else
  if (startsWith(method, "JZS"))            "grey25"  else
  if (startsWith(method, "Zellner"))        "grey60"  else "black"

method_scales <- function(methods, linetype = FALSE) {
  methods <- unique(as.character(methods))
  group <- match(sapply(methods, method_colour_of),
                 c("#0072B2", "#D55E00", "#009E73", "grey25", "grey60", "black"))
  methods <- methods[order(group)]
  list(scale_colour_manual(values = setNames(sapply(methods, method_colour_of), methods),
                           breaks = methods, name = NULL),
       theme(legend.position = "bottom", legend.key.width = unit(2.5, "lines")),
       guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))))
}


## ---- Correlation matrices from the LKJ distribution -------------------------
## Every valid correlation matrix has density proportional to det^(eta - 1);
## eta = 1 is uniform. C-vine sampler of Lewandowski, Kurowicka & Joe (2009):
## draw partial correlations (independent Betas on (-1, 1) whose parameters
## step down by 1/2 per level, ending at Beta(eta, eta)), then convert them
## into ordinary correlations. Used by the motivating simulation
## (R/docs/motivating_simulation.R).
draw_lkj_correlation_matrix <- function(dimension, eta = 1) {
  partial <- matrix(0, dimension, dimension); correlation <- diag(dimension)
  beta_parameter <- eta + (dimension - 1) / 2
  for (level in seq_len(dimension - 1)) {
    beta_parameter <- beta_parameter - 0.5
    for (i in (level + 1):dimension) {
      partial[level, i] <- 2 * rbeta(1, beta_parameter, beta_parameter) - 1
      value <- partial[level, i]
      if (level > 1) for (l in (level - 1):1)
        value <- value * sqrt((1 - partial[l, i]^2) * (1 - partial[l, level]^2)) +
                 partial[l, i] * partial[l, level]
      correlation[level, i] <- correlation[i, level] <- value
    }
  }
  correlation
}


## ---- Whole models: log Bayes factor against the intercept-only model --------
## simple = FALSE returns the log Bayes factor (simple = TRUE overflows to Inf)
log_bayes_factor_model <- function(method, r_squared, n, k) {
  if (k == 0) return(0)
  if (method %in% names(contenders))
    return(log_bayes_factor_vs_intercept_only(r_squared, n, k, prior = contenders[[method]]))
  switch(method,
    `JZS default`  = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = sqrt(2) / 4, simple = FALSE)$bf,
    `Zellner-Siow` = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = 1, simple = FALSE)$bf,
    stop("Unknown method: ", method))
}


## ---- One coefficient --------------------------------------------------------
## Log Bayes factor for one coefficient, every method: the model with it
## against the model without it, each with its own prior. Needs the R-squared
## of both models, n and the number of predictors k in the model with it.
log_bayes_factor_one_coefficient <- function(method, r_squared_full, r_squared_without, n, k)
  log_bayes_factor_model(method, r_squared_full, n, k) -
    log_bayes_factor_model(method, r_squared_without, n, k - 1)

## All methods for one tested column of the predictors
log_bayes_factors_one_coefficient <- function(outcome, predictors, tested_column, methods = model_methods) {
  n <- length(outcome); k <- ncol(predictors)
  r_squared_full    <- sample_r_squared(predictors, outcome)
  r_squared_without <- sample_r_squared(predictors[, -tested_column, drop = FALSE], outcome)
  data.frame(method = methods,
             log_BF10 = vapply(methods, log_bayes_factor_one_coefficient, numeric(1),
                               r_squared_full = r_squared_full, r_squared_without = r_squared_without,
                               n = n, k = k),
             row.names = NULL)
}
