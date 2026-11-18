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
##  Every Bayes factor is a comparison of two models, each with its own
##  prior. For the three priors it is flatcor's: each model fitted with its
##  Stan model, its marginal likelihood estimated by bridge sampling. The
##  references appear in studies 2b and 3b.
## =====================================================================


## ---- Method names and colours (the same in every study) -----------------------------

contenders <- c(
  `flat partial r` = "flat_partial_r",
  `flat R` = "flat_r",
  `flat R-squared` = "flat_r_squared"
)
model_methods <- c(names(contenders), "JZS default", "Zellner-Siow")

## Colour for the figures: one colour per prior, the references in grey
method_colours <- c(
  `flat partial r` = "#0072B2",
  `flat R-squared` = "#009E73",
  `flat R` = "#D55E00",
  `JZS` = "grey25",
  `Zellner` = "grey60"
)
method_colour_of <- function(method) {
  starts <- names(method_colours)[startsWith(method, names(method_colours))]
  if (length(starts) == 0) "black" else method_colours[[starts[1]]]
}

## The colour scale for a figure, with the methods in a fixed order
method_scales <- function(methods) {
  methods <- unique(as.character(methods))
  order_of_colours <- c("#0072B2", "#D55E00", "#009E73", "grey25", "grey60", "black")
  methods <- methods[order(match(sapply(methods, method_colour_of), order_of_colours))]

  list(
    scale_colour_manual(
      values = setNames(sapply(methods, method_colour_of), methods),
      breaks = methods,
      name = NULL
    ),
    theme(legend.position = "bottom", legend.key.width = unit(2.5, "lines")),
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
  )
}


## ---- The Bayes factor between two models -------------------------------------------

## Draws per model for a Bayes factor: ten times the estimation draws (the
## bridgesampling authors' rule), fewer in the smoke test
bayes_factor_draws <- if (identical(Sys.getenv("VALIDATION_SMOKE"), "1")) 1000 else 10000

## The two models as formulas on one data frame, as flatcor takes them: the
## outcome, model A's columns a1, a2, ... and model B's b1, b2, ... (no
## columns: the intercept-only model)
models_as_formulas <- function(predictors_a, predictors_b, outcome) {
  colnames(predictors_a) <- sprintf("a%d", seq_len(ncol(predictors_a)))
  colnames(predictors_b) <- sprintf("b%d", seq_len(ncol(predictors_b)))
  formula_of <- function(columns) {
    reformulate(if (length(columns) == 0) "1" else columns, response = "outcome")
  }
  list(
    formula_a = formula_of(colnames(predictors_a)),
    formula_b = formula_of(colnames(predictors_b)),
    data = data.frame(outcome = outcome, predictors_a, predictors_b)
  )
}

## The log BF of model A (columns predictors_a) against model B
## (predictors_b), for one method, with its approximate Monte Carlo SE and,
## for the three priors, the divergent transitions of the two Stan fits.
## `seed` fixes the Stan fits and bridge sampling (derive it from the
## dataset's ID).
log_bayes_factor_between <- function(method, predictors_a, predictors_b, outcome, seed) {
  # the three priors: flatcor
  if (method %in% names(contenders)) {
    models <- models_as_formulas(predictors_a, predictors_b, outcome)
    result <- suppressWarnings(flatcor::bayes_factor_models(
      models$formula_a,
      models$formula_b,
      models$data,
      prior = contenders[[method]],
      iter_sampling = bayes_factor_draws,
      seed = seed
    ))
    return(c(
      log_bf = result$log_bayes_factor,
      error = result$error,
      divergent = sum(result$divergent)
    ))
  }

  # the references: BayesFactor's regression Bayes factor of each model
  # against the intercept-only model, from R-squared, n and k
  scale <- switch(method,
    `JZS default` = sqrt(2) / 4,
    `Zellner-Siow` = 1,
    stop("Unknown method: ", method)
  )
  log_bf_against_intercept <- function(predictors) {
    if (ncol(predictors) == 0) {
      return(c(bf = 0, error = 0))
    }
    r_squared <- summary(lm(outcome ~ predictors))$r.squared
    bf <- BayesFactor::linearReg.R2stat(
      length(outcome), ncol(predictors), r_squared,
      rscale = scale,
      simple = FALSE
    )
    c(bf = bf$bf, error = bf$properror)
  }
  a <- log_bf_against_intercept(predictors_a)
  b <- log_bf_against_intercept(predictors_b)

  c(
    log_bf = a[["bf"]] - b[["bf"]],
    error = sqrt(a[["error"]]^2 + b[["error"]]^2),
    divergent = 0
  )
}

## A seed for the Stan fits, from a dataset's ID and a label
seed_of <- function(...) abs(digest::digest2int(paste(...))) %% .Machine$integer.max
