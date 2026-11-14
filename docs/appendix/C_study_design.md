# Appendix C. The design of the validation studies (ADEMP)

One description per study of the validation, five studies since the
simplification of 14 November 2026 (see [`REDESIGN.md`](../working/REDESIGN.md) for
the structure and the decisions), in the
ADEMP form of Morris, White and Crowther (2019) as adapted for psychology by
Siepe et al. (2024): **A**ims, **D**ata-generating mechanisms, **E**stimands,
**M**ethods and **P**erformance measures. Every number of simulated datasets
(*n*<sub>sim</sub>) is derived from a target Monte Carlo standard error (MCSE), fixed
here before anything is run. The package under test is flatcor 2.0.0
([dm-k1/flatcor](https://github.com/dm-k1/flatcor), pinned in `renv.lock`).

**One rule for every study:** the three priors are computed by the same
technique and checked by the same routes. flatcor has two implementations,
the same for every prior with only the shape *a* of *R*² ~ Beta(*a*, 1)
changing: the Bayes factor of a model against the intercept-only model (the
likelihood ratio averaged over the prior by `integrate()`; any two models
compare by dividing) and the posterior (one Stan model). A Bayes factor for
one coefficient or one term is the model with it against the model without
it, each with its own prior.

Notation: *n* is the sample size of one simulated dataset, *k* the number of
predictors, *n*<sub>sim</sub> the number of simulated datasets (or prior draws) per
design cell.

---

## How *n*<sub>sim</sub> is chosen

For a performance measure that is a proportion *p* (a rate of decisions,
interval coverage), the Monte Carlo standard error is

$$\mathrm{MCSE} = \sqrt{\frac{p(1-p)}{n_\text{sim}}}, \qquad\text{so}\qquad
n_\text{sim} = \frac{p(1-p)}{\mathrm{MCSE}^2},$$

taken at the worst case *p* = .5 where *p* is not known in advance (Morris et
al., 2019, Table 6). For a mean (a bias, a mean log difference), MCSE =
SD / √*n*<sub>sim</sub>, with the SD taken from the run of 4 November 2026 where one
exists.

| Target | MCSE | *n*<sub>sim</sub> (worst case) |
|---|---|---|
| A proportion reported per design cell | .02 | 625 |
| A proportion reported pooled over cells | .01 | 2,500 in total (met by pooling cells at 625 each) |
| Interval coverage near .95, per cell | .01 | 475 |
| A Kolmogorov–Smirnov check at the .01 level | critical distance ≤ .01 | 1.63² / .01² ≈ 26,600 (effective draws) |

Checks then compare an estimate with its target **relative to its MCSE**
(for example "within 2 MCSE of zero") instead of against a fixed number.

All results are reported with their MCSE, computed by `rsimsum` where it
applies.

---

## Shared data-generating mechanisms

Each dataset has an ID (data set, design cell, replication) and its own
L'Ecuyer-CMRG random stream derived from that ID, so the same ID always gives
the same data, on any machine and in any order. Datasets are regenerated from
their seed, not stored.

| Set | How the data arise | Grid | Used by |
|---|---|---|---|
| **A** | Fixed true effects. Predictors are multivariate normal with unit variances and equal correlation ρ; noise SD 1; slopes in units of the noise SD. | *k* ∈ {2, 5, 10, 20}; *n* ∈ {50, 200, 1,000}; ρ ∈ {0, .3, .7}; other slopes sparse (one at 0.3, the rest 0) or dense (all 0.1); slope of *x*₁ ∈ {0, 0.1, 0.2, 0.3}. That is 72 designs, each at 4 slopes. Extended with a factor of *q* + 1 levels (*q* = 2 to 5; Cohen's *f*² = 0, .02, .15) next to 2 or 10 other predictors that explain *R*² = .5. | 5 |
| **B** | Drawn from each prior: the coefficients and the noise SD from the prior (proper stand-ins for the intercept and noise), then the data. For the Bayes factor's calibration, half the datasets come from the model without *x*₁ (its own prior for *k* − 1 predictors). | *k* ∈ {3, 8}; *n* ∈ {50, 200} (and *n* = 60 at *k* = 8 for calibration of the posterior); one cell per prior. | 3 |
| **D** | *k* growing with *n*: *k* = *n*/4 and *k* = √*n*. Under the null only *R*² is drawn (it is Beta(*k*/2, (*n* − *k* − 1)/2)); with a small true model (3 slopes of 0.3) full datasets are drawn. | *n* ∈ {50, 100, 200, 500, 1,000, 2,000}. | 4 |

---

## The motivating simulation (outside the validation)

`R/docs/motivating_simulation.R`. Not a study: it shows why there are three
priors.

---

## 1. Does the code draw the prior?

`R/studies/01_prior_sampling.R`

- **Aims.** flatcor's Stan model samples *R*, with *R*² ~ Beta(*a*, 1), and a uniform direction. Does it reproduce each prior?
- **Data.** None: the model samples the prior only. Fixed designs with *n* = 60: one predictor, three correlated (.4), eight uncorrelated.
- **Estimands.** The distribution of *r*² (*k* = 1), of *R*² and of each partial correlation.
- **Methods.** The Stan model for each prior.
- **Performance.** The KS distance against the exact distribution (or, for the partial correlations, against the prior as a recipe, `draw_coefficients_from_prior()`), judged against the KS critical value at the .01 level for the effective sample size.
- ***n*<sub>sim</sub>.** Chains lengthened until every checked quantity reaches an ESS of 26,569.
- **Checks, the same four for each prior.** 1.1 *k* = 1: *r*² ~ Beta(*a*, 1); 1.2 *k* = 3: *R*² ~ Beta(*a*, 1); 1.3 *k* = 8: the same; 1.4 *k* = 3: each partial correlation matches the recipe.

## 2. Are the Bayes factors computed correctly?

`R/studies/02_bayes_factors.R`

- **Aims.** Does flatcor's Bayes factor equal the integral it stands for, for every prior?
- **Data.** Grids of *R*², *n* and *k* (*R*² from .001 to .9999, *n* to 10⁵, *k* to 30), and three simulated datasets.
- **Estimands.** The log Bayes factor against the intercept-only model; the log likelihood ratio inside it.
- **Methods.** The package against routes that share no code with it, the same four for every prior: adaptive integration of a second derivation (the prior as a mixture of g-priors); the likelihood ratio the three priors share, against its exact hypergeometric series in 200-bit arithmetic (`Rmpfr`); brute force (draw coefficients from the prior, average how well they predict); bridge sampling on the Stan model (`bridgesampling`).
- **Performance.** The absolute (or, for large evidence, relative) log difference; for the Monte Carlo routes, against their own MCSE.
- ***n*<sub>sim</sub>.** Deterministic, except brute force and bridge sampling, each run until its MCSE is below 0.01.
- **Checks, for each prior.** 2.1 within 1e-5 of `integrate()` (relative 1e-6 where |log BF| > 1); 2.2 the shared likelihood ratio within 1e-6 of its series; 2.3 brute force within 3 MCSE; 2.4 bridge sampling within 3 of its errors.

## 3. Are the posterior and the Bayes factor calibrated?

`R/studies/03_calibration.R`

- **Aims.** Part A: simulation-based calibration of the posterior (Talts et al., 2018), with `SBC`. Part B: the Bayes factor's calibration (Schad et al., 2023): among datasets with a given BF₁₀, the share from the model without *x*₁ equals 1 / (1 + BF₁₀).
- **Data.** Set B. Part A: *k* = 8, *n* = 60 and *k* = 3, *n* = 50. Part B: *k* ∈ {3, 8}, *n* ∈ {50, 200}; generating priors the three and JZS; half the datasets from the model without *x*₁, half from the model with all *k*.
- **Estimands.** Part A: the ranks of the true *R*², *b*₁ and σ. Part B: the observed share without *x*₁ per bin of log BF₁₀ against the predicted share; the logistic calibration slope.
- **Methods.** Part A: the Stan model. Part B: the Bayes factor of the model with *x*₁ against the model without it; on data from a prior that prior's, on JZS data all three.
- **Performance.** Part A: a chi-square over 20 rank bins. Part B: a chi-square over the bins and the logistic slope with its SE.
- ***n*<sub>sim</sub>.** Part A: 1,000 fits per prior and setting. Part B: 20,000 datasets per cell (bins of 625 or more have MCSE ≤ .02).
- **Checks.** 3.1 ranks uniform (*p* > .01), each prior and setting; 3.2 each prior's Bayes factor calibrated on its own prior's data (*p* > .01, slope within 3 SE of 1); 3.3 on JZS data, miscalibration shows (*p* < .001).

## 4. Does each prior meet the objective-Bayes criteria?

`R/studies/04_default_prior_criteria.R`

- **Aims.** The criteria of Bayarri, Berger, Forte and García-Donato (2012) for each prior: predictive matching, information consistency, intrinsic consistency, model-selection consistency (fixed and growing *k*), invariance.
- **Data.** None for the formula checks; simulated *R*² values and datasets for the rates; set D for growing *k*.
- **Estimands.** The slope of the median log Bayes factor for the true model against log *n* (expected: *a* against the intercept-only model, *q*/2 against a model with *q* extra null predictors); for growing *k*, that it is positive.
- **Methods.** The three priors; JZS and Zellner–Siow as references for growing *k*.
- **Performance.** Slopes with bootstrap MCSE, judged by their 95% Monte Carlo interval against a stated tolerance.
- ***n*<sub>sim</sub>.** 2,000 *R*² draws per setting; 400 datasets per setting for growing *k* with a true model.
- **Checks, for each prior.** 4.1 predictive matching (the integral itself at *n* = *k* + 1); 4.2 and 4.3 information consistency; 4.4 to 4.6 the rates; 4.7 consistency under the alternative; 4.8 intrinsic consistency; 4.9 and 4.10 growing *k*; 4.11 invariance under rescaling, shifting and recoding.

## 5. How does each prior behave on data with a known truth?

`R/studies/05_behaviour.R`

- **Aims.** Operating characteristics of the three priors next to JZS and Zellner–Siow. Part A: the Bayes factor for one coefficient (the model with *x*₁ against the model without it): how often it supports either model or stays undecided, how that depends on *k*, how well it separates a null from a small effect. Part B: the same for a factor. Part C: bias, coverage and prediction error of the posterior.
- **Data.** Set A: part A all 72 designs at four slopes; part B the extension with a factor; part C 16 cells (*k* 2 and 10, *n* 50 and 200, ρ 0 and .7, sparse and dense) at slope 0.2, the same datasets as part A's.
- **Estimands.** Parts A and B: the share of each answer at BF 3; the AUC for a null against a small effect and its difference from JZS. Part C: relative bias, coverage of the 95% interval, excess prediction error relative to OLS.
- **Methods.** The three priors (Bayes factors and Stan posteriors), JZS and Zellner–Siow (BayesFactor), and OLS in part C.
- **Performance.** Shares with MCSE; AUC differences and prediction-error ratios with paired bootstrap SEs; bias and coverage with `rsimsum`.
- ***n*<sub>sim</sub>.** Parts A and B: 625 datasets per cell. Part C: 475 per cell.
- **Checks, for every method.** 5.1 the share of true nulls confirmed is the same at *k* = 2 and *k* = 20 (95% Monte Carlo interval of the difference within ±.05); 5.2 JZS for a factor computed precisely enough (proportional error < .05); 5.3 coverage near nominal, each prior; 5.4 prediction at least as good as OLS, each prior; 5.5 no divergent transitions.

---

## Summary

| Study | *n*<sub>sim</sub> per cell | Main cost |
|---|---|---|
| 1 | 26,600 effective draws | Stan, prior only |
| 2 | deterministic; brute force and bridge sampling to MCSE 0.01 | `Rmpfr`, bridge sampling |
| 3 | 1,000 fits (part A); 20,000 datasets (part B) | Stan; the Bayes factors |
| 4 | 2,000 *R*² draws; 400 datasets (growing *k*) | the *k* = 500 regressions |
| 5 | 625 (parts A, B); 475 (part C) | the Bayes factors; BayesFactor; Stan |

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012).
  Criteria for Bayesian model choice with application to variable selection.
  *Annals of Statistics*, *40*(3), 1550–1577.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11),
  2074–2102. https://doi.org/10.1002/sim.8086
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S.
  (2023). Workflow techniques for the robust use of Bayes factors.
  *Psychological Methods*, *28*(6), 1404–1426.
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
  Pawel, S. (2024). Simulation studies for methodological research in
  psychology: A standardized template for planning, preregistration, and
  reporting. *Psychological Methods*. https://doi.org/10.1037/met0000695
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
  Validating Bayesian inference algorithms with simulation-based calibration.
  arXiv:1804.06788.
