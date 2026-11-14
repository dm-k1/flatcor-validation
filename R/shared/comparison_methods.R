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
##  Every Bayes factor is a comparison of two models, each with its own prior.
##  The references appear in study 4 (and later in the behaviour study).
## =====================================================================

## ---- Method names and colours (the same in every study) --------------------
contenders <- c(`flat partial r` = "flat_partial_r", `flat R` = "flat_r",
                `flat R-squared` = "flat_r_squared")
model_methods <- c(names(contenders), "JZS default", "Zellner-Siow")

## Colour for the figures: one colour per prior, the references in grey.
method_colour_of <- function(method)
  if (startsWith(method, "flat partial r")) "#0072B2" else
  if (startsWith(method, "flat R-squared")) "#009E73" else
  if (startsWith(method, "flat R"))         "#D55E00" else
  if (startsWith(method, "JZS"))            "grey25"  else
  if (startsWith(method, "Zellner"))        "grey60"  else "black"

method_scales <- function(methods) {
  methods <- unique(as.character(methods))
  group <- match(sapply(methods, method_colour_of),
                 c("#0072B2", "#D55E00", "#009E73", "grey25", "grey60", "black"))
  methods <- methods[order(group)]
  list(scale_colour_manual(values = setNames(sapply(methods, method_colour_of), methods),
                           breaks = methods, name = NULL),
       theme(legend.position = "bottom", legend.key.width = unit(2.5, "lines")),
       guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))))
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

