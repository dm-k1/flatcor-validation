## =====================================================================
##  Estimates against expectations (the reporting rule, docs/working/REDESIGN.md)
##
##  Every study reports each quantity as one row: what is expected (fixed in
##  advance), the estimate, its Monte Carlo SE and the distance between them
##  in Monte Carlo SEs. No verdict per row: beyond_reference() counts the
##  rows beyond 2 and 3 MCSEs against what chance alone gives if every
##  expectation held (about 4.6% and 0.27%).
## =====================================================================


## Rows of estimates against their expected values; `...` are columns that
## say what each row is (prior, cell, ...)
expectation_rows <- function(..., quantity, expected, estimate, mcse) {
  rows <- data.frame(
    ...,
    quantity = quantity,
    expected = expected,
    estimate = estimate,
    mcse = mcse,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rows$distance_in_mcse <- (rows$estimate - rows$expected) / rows$mcse
  rows
}

## Per group of rows: how many lie beyond 2 and 3 MCSEs of their expected
## value, against what chance alone gives
beyond_reference <- function(rows, by = "quantity") {
  groups <- split(rows, rows[by], drop = TRUE)

  per_group <- lapply(names(groups), function(name) {
    distances <- abs(groups[[name]]$distance_in_mcse)
    data.frame(
      group = name,
      estimates = length(distances),
      beyond_2_mcse = sum(distances > 2, na.rm = TRUE),
      expected_by_chance_2 = round(0.0455 * length(distances), 1),
      beyond_3_mcse = sum(distances > 3, na.rm = TRUE),
      expected_by_chance_3 = round(0.0027 * length(distances), 1)
    )
  })

  do.call(rbind, per_group)
}
