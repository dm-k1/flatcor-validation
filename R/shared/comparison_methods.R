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
## (predictors_b), for one method, with everything computed on the way, as
## one row (a tibble) for the per-dataset files (R/shared/per_dataset.R):
##   id, seed        the name the Stan fits and bridge sampling are seeded
##                   from, and the seed (model B's is seed + 1)
##   log_bf, error   the log BF and its approximate Monte Carlo SE
##   divergent       Stan's divergent transitions of both fits
##   n, k_a, k_b, r_squared_a, r_squared_b
##                   the data's size, both models' numbers of predictors and
##                   their sample R-squared
##   log_ml_a, error_a, ... (and _b)
##                   per model: the log marginal likelihood, its error,
##                   divergences, bridge sampling's iterations, seconds
##   variables       per model, the summary of every Stan variable (flatcor)
##   chains          per model, the sampler per chain (flatcor)
## For JZS and Zellner-Siow: BayesFactor's Bayes factor of each model against
## the intercept-only model and its error, no Stan.
log_bayes_factor_between <- function(method, predictors_a, predictors_b, outcome, name) {
  seed <- seed_of(name)
  sample_r_squared <- function(predictors) {
    if (ncol(predictors) == 0) 0 else summary(lm(outcome ~ predictors))$r.squared
  }
  common <- tibble(
    id = name, seed = seed, method = method, n = length(outcome),
    k_a = ncol(predictors_a), k_b = ncol(predictors_b),
    r_squared_a = sample_r_squared(predictors_a), r_squared_b = sample_r_squared(predictors_b)
  )

  # the three priors: flatcor
  if (method %in% names(contenders)) {
    compile_flatcor_models() # once per worker, one worker at a time (R/shared/pipeline.R)
    models <- models_as_formulas(predictors_a, predictors_b, outcome)
    result <- suppressWarnings(flatcor::bayes_factor_models(
      models$formula_a,
      models$formula_b,
      models$data,
      prior = contenders[[method]],
      iter_sampling = bayes_factor_draws,
      seed = seed
    ))
    per_model <- function(part) {
      m <- result$models[[part]]
      row <- tibble(
        log_ml = m$log_marginal_likelihood, error = m$error, divergent = m$divergent,
        bridge_iterations = m$bridge_sampling$iterations,
        seconds_sampling = m$seconds[["sampling"]],
        seconds_bridge = m$seconds[["bridge_sampling"]]
      )
      names(row) <- paste0(names(row), "_", part)
      row
    }
    long <- function(table) {
      bind_rows(lapply(c("a", "b"), function(part) {
        cbind(model = part, as.data.frame(result$models[[part]][[table]]))
      }))
    }
    return(bind_cols(
      common,
      tibble(log_bf = result$log_bayes_factor, error = result$error,
             divergent = sum(result$divergent)),
      per_model("a"), per_model("b"),
      tibble(variables = list(long("summary")), chains = list(long("sampler")))
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
    r_squared <- sample_r_squared(predictors)
    bf <- BayesFactor::linearReg.R2stat(
      length(outcome), ncol(predictors), r_squared,
      rscale = scale,
      simple = FALSE
    )
    c(bf = bf$bf, error = bf$properror)
  }
  a <- log_bf_against_intercept(predictors_a)
  b <- log_bf_against_intercept(predictors_b)

  bind_cols(common, tibble(
    log_bf = a[["bf"]] - b[["bf"]],
    error = sqrt(a[["error"]]^2 + b[["error"]]^2),
    divergent = 0,
    log_bf_intercept_a = a[["bf"]], error_intercept_a = a[["error"]],
    log_bf_intercept_b = b[["bf"]], error_intercept_b = b[["error"]]
  ))
}
