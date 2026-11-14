## =====================================================================
##  Data set B: data drawn from each prior (docs/appendix/C_study_design.md)
##
##  Used by study 3 (calibration of the posterior and of the Bayes factor).
##  Needs R/data/data_sets.R (with_stream_of()) and flatcor.
##
##  One dataset: predictors with equal correlation .3 (centred, then held
##  fixed: the model conditions on them); the intercept from N(0, 10) and the
##  noise SD from lognormal(0, 1), the proper stand-ins of the Stan model;
##  the coefficients from the generating prior, set up in the predictors'
##  sample covariance (/ n), as the Stan model does. The generating prior is
##  one of the three ("flat_partial_r", "flat_r", "flat_r_squared").
##
##  x1 is drawn with the others, or left out:
##    x1 = "prior"    drawn with the others: data from the model with all k
##                    predictors
##    x1 = "reduced"  0, and x2 ... xk drawn from the generating prior for
##                    those k - 1 predictors (in their own covariance): data
##                    from the model without x1, with its own prior. Half of
##                    the Bayes factor calibration's datasets.
## =====================================================================

sample_covariance <- function(predictors) crossprod(predictors) / nrow(predictors)

## Coefficients from the generating prior, in units of the noise SD
draw_coefficients_B <- function(generating_prior, covariance)
  drop(draw_coefficients_from_prior(1, covariance, prior = generating_prior))

dataset_id_B <- function(generating_prior, k, n, x1, purpose, replication)
  sprintf("B-%s-%s-k%d-n%d-x1_%s-r%d", purpose, generating_prior, k, n, x1, replication)

## One dataset of set B. purpose names the study's use ("sbc",
## "calibration") and is part of the ID, so the uses do not share datasets
## by accident; the same arguments always give the same dataset.
generate_B <- function(generating_prior, k, n, x1 = "prior", purpose, replication,
                       predictor_correlation = 0.3) {
  with_stream_of(dataset_id_B(generating_prior, k, n, x1, purpose, replication), {
    covariance_population <- matrix(predictor_correlation, k, k); diag(covariance_population) <- 1
    predictors <- matrix(rnorm(n * k), n) %*% chol(covariance_population)
    predictors <- scale(predictors, center = TRUE, scale = FALSE)
    attributes(predictors)[c("scaled:center")] <- NULL
    colnames(predictors) <- paste0("x", seq_len(k))
    covariance <- sample_covariance(predictors)
    intercept  <- rnorm(1, 0, 10)
    noise_sd   <- rlnorm(1, 0, 1)
    coefficients <- if (identical(x1, "reduced"))                      # noise units
      c(0, if (k > 1) draw_coefficients_B(generating_prior, covariance[-1, -1, drop = FALSE])) else
      draw_coefficients_B(generating_prior, covariance)
    coefficients <- noise_sd * coefficients                            # data units
    outcome <- intercept + drop(predictors %*% coefficients) + rnorm(n, 0, noise_sd)
    f_squared <- sum((coefficients %*% covariance) * coefficients) / noise_sd^2
    list(predictors = predictors, outcome = outcome, intercept = intercept, noise_sd = noise_sd,
         coefficients = coefficients, r_squared = f_squared / (1 + f_squared))
  })
}

