# What the validation tests

Every study of the validation in one place: its design, what is tested, what is expected and why, with what Monte Carlo error and how many datasets, and the data it uses. It is written before anything is run; the studies report back the same rows, filled in with the estimate and its distance from the expectation in Monte Carlo standard errors (MCSE). No row is judged pass or fail: the table under each result counts how many estimates lie beyond 2 and 3 MCSEs, against what chance alone gives (about 4.6% and 0.3%).

The validation checks three things for each of the three priors (flat partial r, flat R, flat R²), each in two ways:

|   | 1\. Prior sampling | 2\. Posterior estimation | 3\. Bayes factors |
|------------------|------------------|------------------|------------------|
| **a. Computed correctly?** | [1a](#1a-does-the-software-draw-the-prior-it-claims) | [2a](#2a-is-the-posterior-computed-correctly) | [3a](#3a-are-the-bayes-factors-right) |
| **b. Behaves like a default?** | [1b](#1b-is-the-bayes-factor-invariant-to-units-and-coding) | [2b](#2b-how-close-are-the-estimates-to-a-known-truth) | [3b](#3b-do-the-bayes-factors-meet-the-default-prior-criteria) |

Notation: *n* the sample size of one dataset, *k* the number of predictors, *a* the shape of the prior (*R*² \~ Beta(*a*, 1): *k*/2 for flat partial r, ½ for flat R, 1 for flat R²), *n*<sub>sim</sub> the number of datasets per cell. Designs are written as factors crossed: *k* (2, 10) × *n* (50, 200) means the four combinations.

**How many datasets.** For each quantity, *n*<sub>sim</sub> follows from a target Monte Carlo SE fixed in advance (Morris et al., 2019), or, where every dataset costs Stan fits with bridge sampling, from the cost; the column "set by" says which. Numbers set by cost are **provisional**, until the run-time estimate. (Each study also has a batch size: how many datasets form one piece of work, saved as it finishes. That is bookkeeping for resuming and running in parallel, and changes no result.)

------------------------------------------------------------------------

## The data sets

Every dataset is made by one generator in `R/data/`, from its own random stream seeded by its ID (`R/data/building_blocks.R`), so the same ID always gives the same data, on any machine and in any order.

| Data set | File | One dataset | Used by | Why this data |
|------------------|------------------|------------------|------------------|------------------|
| **Prior-only designs** | `prior_only_designs.R` | fixed predictors, *n* = 60: one; three with correlation .4; eight uncorrelated (no outcome: the model draws the prior) | 1a | one predictor (where flat partial r and flat R coincide), a few correlated (the direction is defined relative to their covariance), and more |
| **Example data** | `example_data.R` | *n* = 60, three predictors with correlation .4, slopes 0.3, 0.2, 0, noise SD 1 | 1b | the package's own example: realistic, and the same for every prior |
| **Data from the prior** | `from_the_prior.R` | predictors (correlation .3) × coefficients drawn from the prior being tested + noise; with *x*₁ = "reduced", *x*₁ = 0 and the others from the prior of the model without *x*₁ | 2a, 3a | SBC and Bayes factor calibration are statements about data drawn from the prior itself |
| **Known slopes** | `known_slopes.R` | predictors (correlation ρ) × fixed slopes (*x*₁ = 0.2; the others sparse, one at 0.3, or dense, all 0.1) + noise SD 1 | 2b | estimation is judged against a known truth |
| **Set R²** | `set_r_squared.R` | predictors and an outcome built so that the sample *R*² is exactly a given value | 3b.1, 3b.2 | predictive matching and information consistency are statements about given values of *R*² |
| **Nested models** | `nested_models.R` | standard normal predictors for a smaller and a larger model (the larger adds *q* null predictors, or, for the alternative, a real *x*₁ = 0.15), noise SD 1 | 3b.3–3b.6 | the rate at which the Bayes factor for the true model grows with *n* |
| **Growing *k*** | `growing_k.R` | *k* = *n*/4 or √*n* standard normal predictors (at most 90); no slope (null true), or 3 slopes of 0.3 | 3b.8, 3b.9 | consistency when the number of predictors grows with *n* |

**Data from the prior, in detail.** One dataset = predictors × coefficients + noise:

| Part | How | Shared by the three priors? |
|------------------------|------------------------|------------------------|
| Predictors | *n* × *k*, normal, correlation .3, centred | **yes** |
| Noise | standard normal draws, times the noise SD | **yes** |
| Intercept and noise SD | 2a: drawn from N(0, 10) and lognormal(0, 1), the Stan model's proper stand-ins (σ is a parameter whose rank is checked); 3a: 0 and 1 (the Bayes factor does not depend on them) | **yes** |
| Which model it comes from | 3a: odd replications from the model without *x*₁, even ones from the model with all *k* | **yes** |
| Coefficients | **drawn afresh for every dataset from the prior being tested**: *R*² \~ Beta(*a*, 1) and a uniform direction | **no**: each prior gets effects from itself |

- **No fixed effect size.** SBC and calibration are statements about data drawn from the prior: over all datasets the effect sizes cover everything the prior allows, as often as it allows. With one fixed effect size the check would test something else (a frequentist property, as in 2b).
- **Everything but the coefficients shared** (common random numbers; Morris et al., 2019): only the prior differs between the priors' datasets, so a difference between priors cannot come from different data.
- **Predictor correlation .3.** The checks hold for any predictors, since the model conditions on them. Correlated rather than independent, because each prior's direction is defined relative to the predictors' covariance (the unit conversion, factor 3 of the prior), which independent predictors would barely exercise; .3 is a modest, realistic value (0 or .5 would serve as well).

------------------------------------------------------------------------

## 1a. Does the software draw the prior it claims?

flatcor's Stan model with the data switched off should draw exactly the prior.

**Design:** prior (3) × predictor design (one, three correlated, eight) = 9 settings. Data: the prior-only designs.

| Quantity | Expected | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|
| 5th, 25th, 50th, 75th, 95th percentile of *r*² (one predictor) | those of Beta(*a*, 1) | the prior is *R*² \~ Beta(*a*, 1) | `posterior::mcse_quantile()` (allows for autocorrelated draws) | 26,569 effective draws per setting | target: a KS distance of .01 detectable at the .01 level |
| the same of *R*² (three, eight predictors) | those of Beta(*a*, 1) | the same | the same | the same | the same |
| the same of each partial correlation (three predictors) | those of the prior as a recipe (`draw_coefficients_from_prior()`) | a second, independent route to the same prior | the same, plus the recipe's own | the same | the same |

Also shown: the KS distance of each quantity to its expected distribution, next to its .01 critical value, as a description of the whole distribution.

------------------------------------------------------------------------

## 1b. Is the Bayes factor invariant to units and coding?

Bayarri et al. (2012): rescaling, shifting or recoding the predictors and the outcome leaves the Bayes factor unchanged.

**Design:** prior (3) × change (rescaled: *x*₁ → 100 *x*₁ + 7, *x*₂ → −0.01 *x*₂, *x*₃ → *x*₃ − 50, *y* → 1000 *y* − 3; recoded: *x*₁, *x*₂ → *x*₁ + *x*₂, *x*₁ − *x*₂) × Bayes factor (whole model against intercept only; *x*₁ with against without). Data: the example data.

| Quantity | Expected | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|
| **pooled**: the mean change in log BF over the datasets, per prior, change and Bayes factor | 0 | the prior is built in units of explained variance | √(sum of the datasets' MCSE²) / *n*<sub>sim</sub> | 20 datasets per prior | target: with an error of about 0.01 per change, the mean has an MCSE of about 0.003, so a systematic change of 1% in the Bayes factor stands out at about 3 MCSE |
| per dataset: the change in log BF | 0 | the same | bridge sampling's approximate errors of both Bayes factors | the same | the same |

------------------------------------------------------------------------

## 2a. Is the posterior computed correctly?

Simulation-based calibration (Talts et al., 2018): draw the parameters from the prior, data from them, fit; if the posterior is right, the rank of each true value among its posterior draws is uniform.

**Design:** prior (3) × setting (*k* = 8, *n* = 60; *k* = 3, *n* = 50) = 6 cells. Data: data from the prior, intercept and noise SD drawn.

| Quantity, for *R*², *b*₁ and σ | Expected | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|
| mean of the scaled rank (0–1) | ½ | uniform ranks; a bias shifts them | SD of a uniform rank / √*n*<sub>sim</sub> | 1,000 fits per cell | target: a miscalibration of about 5% of the rank distribution |
| SD of the scaled rank | that of a uniform rank, about .289 | a posterior too wide or too narrow changes it | the SD of a sample SD of uniform ranks | the same | the same |
| chi-square over 20 rank bins | 19 | uniform ranks | its SD under uniformity, √38 | the same | the same |

Also shown: Stan's divergent transitions (counts) and the SBC package's ECDF plots.

------------------------------------------------------------------------

## 2b. How close are the estimates to a known truth?

How each prior's posterior estimates behave on data with fixed slopes. The references are not guaranteed for a Bayesian estimate, so this study describes.

**Design:** *k* (2, 10) × *n* (50, 200) × predictor correlation (0, .7) × other slopes (sparse, dense) = 16 cells; slope of *x*₁ = 0.2. Methods: the three priors (Stan), JZS and Zellner–Siow (BayesFactor), OLS. Data: known slopes.

| Quantity, for *x*₁ | Reference | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|
| coverage of the 95% interval | .95 | nominal | √(*p*(1 − *p*)/*n*<sub>sim</sub>) | 475 datasets per cell | target: a coverage near .95 with MCSE .01 needs .95 × .05 / .01² = 475 |
| bias of the posterior mean, as a share of 0.2 | 0 | OLS is unbiased; a prior that shrinks pulls it below 0 | SD / √*n*<sub>sim</sub> | the same | the same |
| excess prediction error / OLS's | 1 | below 1 means better than OLS | paired bootstrap | the same | the same |

------------------------------------------------------------------------

## 3a. Are the Bayes factors right?

**Design:** prior (3) × setting (*k* = 3, *n* = 50; *k* = 8, *n* = 200) = 6 cells; each prior generates its own data and computes its own Bayes factor for *x*₁. Data: data from the prior, intercept 0 and noise SD 1, half the datasets from the model without *x*₁.

| Part | Quantity | Expected | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|------------|
| A | share of datasets without *x*₁, per bin of log BF₁₀ | the mean of 1 / (1 + BF₁₀) in the bin | **calibration** (Schad et al., 2023): among datasets drawn from the prior, P(without *x*₁ \| BF₁₀) = 1 / (1 + BF₁₀), if the Bayes factors are right | the Monte Carlo SD of the share under calibration | 1,000 datasets per cell | cost (**provisional**): about .05 per bin; .02 per bin would need about 10,000 |
| A | calibration intercept | 0 | the same, as a logistic line of "without *x*₁" on log BF₀₁ | SE of the coefficient | the same | cost: SE about .07 |
| A | calibration slope | 1 | the same | SE of the coefficient | the same | cost: SE about .07 (a slope of .85 at about 2 SE) |
| B | spread of repeated log BFs / reported error | 1 | bridgesampling's approximate error should be the real spread (Gronau et al., 2020: repeat with fresh posterior draws) | relative SE of a pooled SD | 20 datasets × 5 runs per cell | target: relative SE about .08 (a ratio of 1.4 at about 5 SE) |

Also shown: bridge sampling's approximate errors and Stan's divergent transitions. Bins of log BF₁₀ run from −2.5 to 6, finer near 0, where most datasets fall; bins with fewer than 10 datasets are left out.

------------------------------------------------------------------------

## 3b. Do the Bayes factors meet the default-prior criteria?

The criteria of Bayarri, Berger, Forte and García-Donato (2012); invariance is [1b](#1b-is-the-bayes-factor-invariant-to-units-and-coding). Where the expectation is a direction or a limit, the row says so.

**Designs:**

- 3b.1: prior (3) × *k* (1, 2, 3) at *n* = *k* + 1 = 9 cells. Data: set R².
- 3b.2: prior (3) × *R*² (.9, .99, .999, .9999) at *n* = 20, *k* = 3 = 12 cells. Data: set R².
- 3b.3–3b.6: setting (3b.3: one null next to two real; 3b.4: *x*₁ real, *q* = 1, 3, 6; 3b.5: intercept only, *q* = 1, 3, 6; 3b.6: alternative) × *n* (50, 200, 1,000, 5,000; the alternative from 200) × prior (3) = 93 cells. Data: nested models.
- 3b.7: prior (3) × *n* (30, 300, 100,000) = 9 cells; the prior as a recipe, no Bayes factors.
- 3b.8, 3b.9: growth (*k* = *n*/4 at *n* = 50, 100, 200, 360; *k* = √*n* at *n* = 50 to 2,000) × truth (null; a small model) × method (the three priors, JZS, Zellner–Siow) = 100 cells. Data: growing *k*.

| Item | Quantity | Expected | Why | MCSE | *n*<sub>sim</sub> | Set by |
|------------|------------|------------|------------|------------|------------|------------|
| 3b.1 | log BF against the intercept-only model at *n* = *k* + 1 | 0 | predictive matching: every model fits exactly, the data say nothing | bridge sampling's error | one dataset per cell | the property is about one given dataset |
| 3b.2 | log BF as *R*² → 1 | increases (a direction) | information consistency | the same | the same | the same |
| 3b.3 | slope of the median log BF for the truth against log *n* | ½ (asymptotic) | model-selection consistency: evidence for a true null grows like *q*/2 log *n* | bootstrap | 100 datasets per cell | cost (**provisional**): slope MCSE about .05 to .1 |
| 3b.4 | the same | *q*/2 (asymptotic) | the same | the same | the same | the same |
| 3b.5 | the same, against the intercept-only model | *a* (asymptotic) | the rate against the intercept-only model is the prior's shape | the same | the same | the same |
| 3b.6 | median log BF₁₀ for a real effect | increases with *n* (a direction) | consistency under the alternative | the same | the same | the same |
| 3b.7 | KS distance of the prior with the sample covariance to that with the population covariance | falls below its .01 critical value (a limit) | intrinsic consistency | — | 26,569 draws | target: as in 1a |
| 3b.8, 3b.9 | slope of the median log BF for the truth against log *n* | positive (a direction) | consistency as *k* grows | bootstrap | 50 datasets per cell | cost (**provisional**): enough to see the sign |

Also shown: Stan's divergent transitions.

**Choices.** *n* from 50 to 5,000 for the rates: a slope against log *n* needs a wide range. *k* at most 90 for growing *k*: a model has *k* + 3 parameters, and bridge sampling has been shown to work up to 96 (Gronau et al., 2017). Not included: whether information consistency starts exactly at *n* = *k* + 3 (it needs *R*² within 10⁻⁸ of 1, beyond what a sampler resolves).

------------------------------------------------------------------------

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria for Bayesian model choice with application to variable selection. *Annals of Statistics*, *40*(3), 1550–1577.
- Gronau, Q. F., et al. (2017). A tutorial on bridge sampling. *Journal of Mathematical Psychology*, *81*, 80–97.
- Gronau, Q. F., Singmann, H., & Wagenmakers, E.-J. (2020). bridgesampling: An R package for estimating normalizing constants. *Journal of Statistical Software*, *92*(10), 1–29.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11), 2074–2102.
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S. (2023). Workflow techniques for the robust use of Bayes factors. *Psychological Methods*, *28*(6), 1404–1426.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788.
