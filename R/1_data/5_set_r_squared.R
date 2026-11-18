## =====================================================================
##  Data with a set sample R-squared (study 3b, items 3b.1 and 3b.2). Why:
##  docs/validation_tests.md, "The data sets".
##
##  Predictors drawn; the outcome built from a fitted part and a residual
##  orthogonal to it, so that its sample R-squared is exactly `r_squared`
##  (at n = k + 1 there is no residual: R-squared is 1).
## =====================================================================


data_set_r_squared <- function(n, k, r_squared) {
  with_stream_of(sprintf("set-r-squared-k%d-n%d", k, n), {
    x <- correlated_predictors(n, k, 0)
    centred <- sweep(x, 2, colMeans(x))

    fitted <- drop(centred %*% rnorm(k))
    fitted <- fitted / sqrt(sum(fitted^2))

    residual <- if (n > k + 1) qr.resid(qr(cbind(1, centred)), rnorm(n)) else rep(0, n)
    if (sum(residual^2) > 0) residual <- residual / sqrt(sum(residual^2))

    list(x = x, y = sqrt(r_squared) * fitted + sqrt(1 - r_squared) * residual)
  })
}
