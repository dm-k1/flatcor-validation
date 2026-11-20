# Per-dataset results (layer 1)

One row per dataset (and prior or method) of every study, with everything
the run computed: the numbers the metrics in `../studies/` are made from,
and what else came out on the way. Written by the pipeline, one step per
file (`R/shared/per_dataset.R`); the metrics are then computed from these
files as they are written here, so they can be redone from them.

The validation's metrics (layer 2) are exactly those of
`docs/validation_tests.md`. Anything else (comparisons between priors,
other summaries) is a secondary analysis of these files, outside the
validation.

**Columns every Bayes factor row has** (1b, 3a, 3b; `R/shared/comparison_methods.R`):
`id` (the name its seed comes from: the dataset can be made again from it),
`seed` (model B's is `seed + 1`), `n`, `k_a`, `k_b`, `r_squared_a`,
`r_squared_b` (the data's sample R² under each model), `log_bf`, `error`
(bridge sampling's approximate SE of the log BF), `divergent`; per model
(`_a`, `_b`): `log_ml`, `error`, `divergent`, `bridge_iterations`,
`seconds_sampling`, `seconds_bridge`; then `seconds` (the whole row) and
`machine`. JZS and Zellner–Siow rows have BayesFactor's Bayes factor of each
model against the intercept-only model instead (`log_bf_intercept_a`, ...).

**Files `_variables.csv.gz` and `_chains.csv.gz`** next to a study's file
hold, joined by `id` and `model` (`a`, `b`): the summary of every Stan
variable (mean, median, sd, mad, quantiles, rhat, ess_bulk, ess_tail) and
the sampler per chain (divergences, treedepth hits, E-BFMI, step size,
leapfrog steps, warmup and sampling seconds), as flatcor returns them.

| File | One row is |
|---|---|
| `study_1a_prior_draws.csv.gz` | one prior draw of a quantity: `prior`, `design`, `quantity`, `chain`, `iteration`, `value` |
| `study_1a_recipe_draws.csv.gz` | one draw of the prior as a recipe (partial correlations, three-predictor design) |
| `study_1a_settings.csv` (+ `_chains`) | one setting's run: seed, iterations per chain, attempts, smallest ESS, seconds |
| `study_1b_log_bayes_factors.csv` (+ both) | one of a dataset's five log BFs: `bayes_factor` (whole model / x1), `version` (original, rescaled, recoded) |
| `study_2a_sbc_fits.csv` | one variable of one SBC fit: true value, rank, posterior summary, rhat, ESS (the SBC package) |
| `study_2a_sbc_diagnostics.csv` | one SBC fit's sampler diagnostics (the SBC package) |
| `study_2b_fits.csv` (+ both) | one method's fit of one dataset: x1's estimate and interval, prediction error; `_variables` has every parameter (OLS: its coefficient table with t and p) |
| `study_3a_calibration.csv` (+ both) | one dataset's log BF10, with the true values it was drawn with |
| `study_3a_repeated_runs.csv` (+ both) | one of several runs on the same dataset (`run`) |
| `study_3b_formula.csv` (+ both) | one prior at one small *n* (3b.1, 3b.2) |
| `study_3b_rates.csv` (+ both) | one dataset's log BF (3b.3–3b.6): `log_bf` for the truth, `log_bf_smaller` of the smaller model against the larger |
| `study_3b_growing_k.csv` (+ both, one per growth rule: `_k_n_4`, `_k_sqrt_n`) | one method on one dataset as *k* grows with *n* (3b.8, 3b.9); its `_variables` and `_chains` are split by `growth` to stay under GitHub's 100 MB per file |
| `study_3b_intrinsic_draws.csv.gz` | one prior draw of x1's coefficient with the sample or the population covariance (3b.7) |

Numbers are written with 15 significant digits (`data.table::fwrite`).
