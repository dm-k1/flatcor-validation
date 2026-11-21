# What the validation tests

Every study of the validation in one place: its design, what is tested, what is expected and why, with what Monte Carlo error and how many datasets, and the data it uses. It is written before anything is run; the studies report back the same rows, filled in with the estimate and its distance from the expectation in Monte Carlo standard errors (MCSE). No row is judged pass or fail: the table under each result counts how many estimates lie beyond 2 and 3 MCSEs, against what chance alone gives (about 4.6% and 0.3%).

The validation checks three things for each of the three priors (flat partial r, flat R, flat R²), each in two ways:

|   | 1\. Prior sampling | 2\. Posterior estimation | 3\. Bayes factors |
|------------------|------------------|------------------|------------------|
| **a. Computed correctly?** | [1a](#1a-does-the-software-draw-the-prior-it-claims) | [2a](#2a-is-the-posterior-computed-correctly) | [3a](#3a-are-the-bayes-factors-right) |
| **b. Behaves like a default?** | [1b](#1b-is-the-bayes-factor-invariant-to-units-and-coding) | [2b](#2b-how-close-are-the-estimates-to-a-known-truth) | [3b](#3b-do-the-bayes-factors-meet-the-default-prior-criteria) |

Notation: *n* the sample size of one dataset, *k* the number of predictors, *a* the shape of the prior (*R*² \~ Beta(*a*, 1): *k*/2 for flat partial r, ½ for flat R, 1 for flat R²), *n*<sub>sim</sub> the number of datasets per cell. Designs are written as factors crossed: *k* (2, 10) × *n* (50, 200) means the four combinations.

**How many datasets.** By the standard procedure for simulation studies (Burton et al., 2006; Morris et al., 2019), every quantity gets an accuracy Δ in advance: the smallest deviation from its expectation that should show. Its Monte Carlo SE must be at most Δ / 1.96, so that a deviation of Δ reaches the 95% line; equivalently *n*<sub>sim</sub> = (1.96 σ / Δ)², with σ the SD of the quantity's estimate over datasets. (Burton et al.'s optional power term is left out: 95%, without power.)

Δ follows the field's conventions, at their looser end throughout, so that the whole validation fits a budget of about €20 of computing:

- **coverage**: a Monte Carlo SE of 1% (Morris et al., 2019, use 0.4–0.5% in their examples; 1% is the common upper end), so Δ = .0196;
- **every other quantity with an expected value**: a relative error of 10% of that value (the common criteria for relative bias are 5% and 10%); for a log Bayes factor, a Bayes factor off by 10% (log 1.1 = .095); for an expectation of 0, the same scale as its companion quantity, said in its row; for the percentiles of a distribution (1a), 10% of that distribution's SD, its natural unit (10% of a percentile itself would demand the impossible for a median of 0, as the partial correlations have);
- **a direction** ("increases", "positive"), whose size is not known in advance: the same 10%, relative to the slope itself, with a floor. Its Monte Carlo SE must be at most max(.1 × |estimate|, Δ) / 1.96, with Δ = .05 per log *n* (10% of the smallest rate in the same study, ½). This is the relative fixed-width rule (Glynn & Whitt, 1992; for MCMC output, Flegal & Gong, 2015): a slope of 30 is then known to within 10%, as every other quantity, instead of to within .05. The floor ends the rule for a slope near 0, where 10% of the estimate would demand the impossible. The rule asks for the same precision whether a slope is positive or negative, so it does not stop on a slope that looks clearly positive (which would be optional stopping).

Where σ is known in advance (coverage, the ranks of 2a), *n*<sub>sim</sub> is computed. Elsewhere it is **measured**, by the classical sequential fixed-width rule (Chow & Robbins, 1965; shown valid for stochastic simulations by Glynn & Whitt, 1992, and for MCMC output by Jones et al., 2006): every cell starts with a **first round of 50 datasets**, fixed in advance, so that its Monte Carlo SE is estimated from enough datasets (the rule's guard against stopping on a lucky, too small estimate; 50 is our choice, as bootstrap SEs of medians and slopes are unreliable on fewer); after each round the SE is estimated again, and the cell **stops at the first round whose Monte Carlo SE is at most Δ / 1.96**. The rule looks only at the Monte Carlo SE, never at the estimate itself, and it is fixed before the run: adapting a simulation's size becomes a questionable practice when it is undisclosed or serves a preferred result (Pawel et al., 2024). The whole run has a budget of 2,500 core hours (about €20; `docs/working/run_time_estimate.md`). A round that would take the run past it is not started: the run stops before it and reports what it would cost (from what a dataset of each study has cost so far). The researcher then decides, on that cost alone: raise the budget and continue, or stop, and report every cell that has not reached its target as such, with its Monte Carlo SE. This is written down before the run, as the template of Siepe et al. (2024) asks.

**Monte Carlo errors are measured, never taken from a method's own report** (Koehler et al., 2009; Morris et al., 2019). A quantity averaged over datasets gets the SD of its estimates over those datasets / √*n*<sub>sim</sub> (or a bootstrap, for medians and slopes), which includes every source of noise, bridge sampling's among them. Only a quantity about one given dataset (1b per dataset, 3b.1, 3b.2) needs a method's error: there, bridge sampling's approximate error is multiplied by the factor by which it falls short of the real spread, measured in 3a part B (Gronau et al., 2020, recommend measuring it from repeated runs with fresh posterior draws). 3a part B therefore runs before these quantities are judged.

(Each study also has a batch size: how many datasets form one piece of work, saved as it finishes. That is bookkeeping for resuming and running in parallel, and changes no result.)

------------------------------------------------------------------------

## The data sets

Every dataset is made by one generator in `R/1_data/` (numbered in the order the studies first use them), from its own random stream seeded by its ID (`R/shared/random_streams.R`), so the same ID always gives the same data, on any machine and in any order.

**Randomness.** Every random number in the validation follows from one seed, `root_seed` in `R/shared/random_streams.R` (it is 0), and a name. Everything random has a name that says what it is: a dataset's name is its data set, cell and replication; a Stan fit's name is its dataset's name and the prior; a bootstrap's name is what it bootstraps. The name and the root seed are hashed (`digest::digest2int`) to a seed, used in one of two ways:

- R's random numbers (making data, bootstraps, BayesFactor's posterior draws) are drawn inside `with_stream_of(name, ...)`: R gets the seed, draws, and is then put back as it was, so making one dataset never changes another. A bootstrap that gives a cell's Monte Carlo SE has a stream of its own per cell, named after the cell, so a cell's SE depends only on its own datasets, not on how many the other cells have.
- Stan and bridge sampling get the seed as `seed = seed_of(...)`. flatcor uses it for Stan (the second model of a Bayes factor gets seed + 1) and for bridge sampling's own random numbers.

targets' own seed, with which it seeds R before each step, is the root seed too, so nothing random is left outside it. The same name always gives the same numbers, on any machine, in any order, on any number of workers and in any round. A new round adds new replication numbers, so new names: finished work is never drawn again, and a piece of work computed again gives the same result. The smoke test and the dry run use the same names for the corner of the design they run. **To run the validation again with fresh random numbers,** change `root_seed` and commit: every step then reruns, and every results file records the root seed it was made with. With 0, every number is the same as in v1 (the hash's seed was 0 there). Two details: in 2a the fits of one batch share one Stan seed (SBC takes one seed per batch), so changing 2a's batch size changes its seeds; and the run-time estimate's reference task (`R/timing/`) keeps its own fixed seed, as a yardstick that must be the same on every machine. Results can still differ in the last digits between processors with the same seeds.

| Data set | File | One dataset | Used by | Why this data |
|------------------|------------------|------------------|------------------|------------------|
| **Prior-only designs** | `1_prior_only_designs.R` | fixed predictors, *n* = 60: one; three with correlation .4; eight uncorrelated (no outcome: the model draws the prior) | 1a | one predictor (where flat partial r and flat R coincide), a few correlated (the direction is defined relative to their covariance), and more |
| **Ordinary regression** | `2_ordinary_regression.R` | *n* = 60, three predictors with correlation .4, slopes 0.3, 0.2, 0, noise SD 1 | 1b | the package's own example: realistic, and the same for every prior |
| **Data from the prior** | `3_from_the_prior.R` | predictors (correlation .3) × coefficients drawn from the prior being tested + noise; with *x*₁ = "reduced", *x*₁ = 0 and the others from the prior of the model without *x*₁ | 2a, 3a | SBC and Bayes factor calibration are statements about data drawn from the prior itself |
| **Known slopes** | `4_known_slopes.R` | predictors (correlation ρ) × fixed slopes (*x*₁ = 0.2; the others sparse, one at 0.3, or dense, all 0.1) + noise SD 1 | 2b | estimation is judged against a known truth |
| **Set R²** | `5_set_r_squared.R` | predictors and an outcome built so that the sample *R*² is exactly a given value | 3b.1, 3b.2 | predictive matching and information consistency are statements about given values of *R*² |
| **Nested models** | `6_nested_models.R` | standard normal predictors for a smaller and a larger model (the larger adds *q* null predictors, or, for the alternative, a real *x*₁ = 0.15), noise SD 1 | 3b.3–3b.6 | the rate at which the Bayes factor for the true model grows with *n* |
| **Growing *k*** | `7_growing_k.R` | *k* = *n*/4 or √*n* standard normal predictors (at most 90); no slope (null true), or 3 slopes of 0.3 | 3b.8, 3b.9 | consistency when the number of predictors grows with *n* |

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

| Quantity | Expected | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|
| 5th, 25th, 50th, 75th, 95th percentile of *r*² (one predictor) | those of Beta(*a*, 1) | the prior is *R*² \~ Beta(*a*, 1) | `posterior::mcse_quantile()` (allows for autocorrelated draws) | 10% of the quantity's SD under the prior | measured: the chains are lengthened until every percentile reaches its target |
| the same of *R*² (three, eight predictors) | those of Beta(*a*, 1) | the same | the same | the same | the same |
| the same of each partial correlation (three predictors) | those of the prior as a recipe (`draw_coefficients_from_prior()`) | a second, independent route to the same prior | the same, plus the recipe's own | the same | the same |

Also shown: the KS distance of each quantity to its expected distribution, next to its .01 critical value, as a description of the whole distribution.

------------------------------------------------------------------------

## 1b. Is the Bayes factor invariant to units and coding?

Bayarri et al. (2012): rescaling, shifting or recoding the predictors and the outcome leaves the Bayes factor unchanged.

**Design:** prior (3) × change (rescaled: *x*₁ → 100 *x*₁ + 7, *x*₂ → −0.01 *x*₂, *x*₃ → *x*₃ − 50, *y* → 1000 *y* − 3; recoded: *x*₁, *x*₂ → *x*₁ + *x*₂, *x*₁ − *x*₂) × Bayes factor (whole model against intercept only; *x*₁ with against without). Data: the ordinary regression data.

| Quantity | Expected | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|
| **pooled**: the mean change in log BF over the datasets, per prior, change and Bayes factor | 0 | the prior is built in units of explained variance | SD of the changes over the datasets / √*n*<sub>sim</sub> | .095: the Bayes factor off by 10% (log 1.1) | measured (datasets per prior) |
| per dataset: the change in log BF | 0 | the same | bridge sampling's approximate errors of both Bayes factors, each × the shortfall factor of 3a part B | — (one dataset; descriptive) | the same |

------------------------------------------------------------------------

## 2a. Is the posterior computed correctly?

Simulation-based calibration (Talts et al., 2018): draw the parameters from the prior, data from them, fit; if the posterior is right, the rank of each true value among its posterior draws is uniform.

**Design:** prior (3) × setting (*k* = 8, *n* = 60; *k* = 3, *n* = 50) = 6 cells. Data: data from the prior, intercept and noise SD drawn.

| Quantity, for *R*², *b*₁ and σ | Expected | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|
| mean of the scaled rank (0–1) | ½ | uniform ranks; a bias shifts them | SD of a uniform rank / √*n*<sub>sim</sub> | .05: 10% of ½ | 131 fits per cell: (.2915 × 1.96 / .05)², .2915 being the SD of a rank from 0 to 100 (ranks thinned to 100 draws) |
| SD of the scaled rank | that of a uniform rank, about .289 | a posterior too wide or too narrow changes it | the SD of a sample SD of uniform ranks | .029: 10% of .289 | 76 needed; 131 cover it |
| chi-square over 20 rank bins | 19 | uniform ranks | its SD under uniformity, √38 | — (an overall check) | the same 131 |

Also shown: Stan's divergent transitions (counts) and the SBC package's ECDF plots.

------------------------------------------------------------------------

## 2b. How close are the estimates to a known truth?

How each prior's posterior estimates behave on data with fixed slopes. The references are not guaranteed for a Bayesian estimate, so this study describes.

**Design:** *k* (2, 10) × *n* (50, 200) × predictor correlation (0, .7) × other slopes (sparse, dense) = 16 cells; slope of *x*₁ = 0.2. Methods: the three priors (Stan), JZS and Zellner–Siow (BayesFactor), OLS. Data: known slopes.

| Quantity, for *x*₁ | Reference | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|
| coverage of the 95% interval | .95 | nominal | √(*p*(1 − *p*)/*n*<sub>sim</sub>) | .0196: an MCSE of 1%, the coverage convention | 475 per cell: .95 × .05 / .01² |
| bias of the posterior mean, as a share of 0.2 | 0 | OLS is unbiased; a prior that shrinks pulls it below 0 | SD / √*n*<sub>sim</sub> | .02: 10% of the slope 0.2 | measured |
| excess prediction error / OLS's | 1 | below 1 means better than OLS | paired bootstrap | .1: 10% of OLS's | measured |

------------------------------------------------------------------------

## 3a. Are the Bayes factors right?

**Design:** prior (3) × setting (*k* = 3, *n* = 50; *k* = 8, *n* = 200) = 6 cells; each prior generates its own data and computes its own Bayes factor for *x*₁. Data: data from the prior, intercept 0 and noise SD 1, half the datasets from the model without *x*₁.

| Part | Quantity | Expected | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|------------|
| A | **data-averaged posterior**: the mean over the datasets of P(without *x*₁ \| data) = 1 / (1 + BF₁₀) | ½, the prior probability of the model without *x*₁ | **calibration** (Schad et al., 2023): if the Bayes factors are right, the posterior model probabilities, averaged over data drawn from the prior, return the prior model probabilities | SD over the datasets / √*n*<sub>sim</sub> | .05: 10% of ½ | measured |
| A | share of datasets without *x*₁, per bin of log BF₁₀ | the mean of 1 / (1 + BF₁₀) in the bin | the same, per Bayes factor: P(without *x*₁ \| BF₁₀) = 1 / (1 + BF₁₀) | the Monte Carlo SD of the share under calibration | — (descriptive; the figure) | those of the data-averaged posterior |
| A | calibration intercept and slope | 0 and 1 | the same, as a logistic line of "without *x*₁" on log BF₀₁ | SE of the coefficient | — (descriptive) | the same |
| B | spread of repeated log BFs / reported error | 1 | bridgesampling's approximate error should be the real spread (Gronau et al., 2020: repeat with fresh posterior draws) | relative SE of a pooled SD | .1: 10% of 1 | measured (datasets of 5 runs each) |

Also shown: bridge sampling's approximate errors and Stan's divergent transitions. Bins of log BF₁₀ run from −2.5 to 6, finer near 0, where most datasets fall; bins with fewer than 10 datasets are left out.

**Why the data-averaged posterior decides, and the line describes.** Both test the same thing, but they learn from different datasets. The data-averaged posterior uses every dataset's posterior probability in full, also an extreme one (.001 or .999). A logistic line learns only from datasets whose Bayes factor is near 1, as only there is it uncertain which model made the data; where Bayes factors are extreme, few datasets inform it. A first round of 50 datasets per cell gave the line SEs of .38–1.35: a Δ of .1 for its intercept and slope would ask for 2,800 to 35,000 datasets per cell (63,000 in all; flat partial r at *k* = 8 alone 35,000, as only 4 of its 50 datasets had |log BF₁₀| < 1.5), against 130 to 300 per cell for the data-averaged posterior. The data-averaged posterior cannot see Bayes factors too loud on one side and too timid on the other, which cancel in the mean; the bins and the line show that pattern, and are reported with their Monte Carlo errors. Repeating a dataset's computation does not help the line: bridge sampling's noise (about .005 on the log Bayes factor) is far below the line's, which comes from each dataset adding a single yes or no.

------------------------------------------------------------------------

## 3b. Do the Bayes factors meet the default-prior criteria?

The criteria of Bayarri, Berger, Forte and García-Donato (2012); invariance is [1b](#1b-is-the-bayes-factor-invariant-to-units-and-coding). Where the expectation is a direction or a limit, the row says so.

**Designs:**

- 3b.1: prior (3) × *k* (1, 2, 3) at *n* = *k* + 1 = 9 cells. Data: set R².
- 3b.2: prior (3) × *R*² (.9, .99, .999, .9999) at *n* = 20, *k* = 3 = 12 cells. Data: set R².
- 3b.3–3b.6: setting (3b.3: one null next to two real; 3b.4: *x*₁ real, *q* = 1, 3, 6; 3b.5: intercept only, *q* = 1, 3, 6; 3b.6: alternative) × *n* (50, 200, 1,000, 5,000; the alternative from 200) × prior (3) = 93 cells. Data: nested models.
- 3b.7: prior (3) × *n* (30, 300, 100,000) = 9 cells; the prior as a recipe, no Bayes factors.
- 3b.8, 3b.9: growth (*k* = *n*/4 at *n* = 50, 100, 200, 360; *k* = √*n* at *n* = 50 to 2,000) × truth (null; a small model) × method (the three priors, JZS, Zellner–Siow) = 100 cells. Data: growing *k*.

| Item | Quantity | Expected | Why | MCSE | Δ | *n*<sub>sim</sub> |
|------------|------------|------------|------------|------------|------------|------------|
| 3b.1 | log BF against the intercept-only model at *n* = *k* + 1 | 0 | predictive matching: every model fits exactly, the data say nothing | bridge sampling's error × the shortfall factor of 3a part B | — (about one given dataset) | one per cell |
| 3b.2 | log BF as *R*² → 1 | increases (a direction) | information consistency | the same | the same | the same |
| 3b.3 | slope of the median log BF for the truth against log *n* | ½ (asymptotic) | model-selection consistency: evidence for a true null grows like *q*/2 log *n* | bootstrap | 10% of ½: .05 | measured |
| 3b.4 | the same | *q*/2 (asymptotic) | the same | the same | 10% of *q*/2 | measured |
| 3b.5 | the same, against the intercept-only model | *a* (asymptotic) | the rate against the intercept-only model is the prior's shape | the same | 10% of *a* | measured |
| 3b.6 | median log BF₁₀ for a real effect | increases with *n* (a direction) | consistency under the alternative | the same | 10% of the slope, at least .05 per log *n* (a direction: relative rule) | measured |
| 3b.7 | KS distance of the prior with the sample covariance to that with the population covariance | falls below its .01 critical value (a limit) | intrinsic consistency | — | — (a limit) | 26,569 draws: a KS distance of .01 detectable at the .01 level (draws only, no Stan fits) |
| 3b.8, 3b.9 | slope of the median log BF for the truth against log *n* | positive (a direction) | consistency as *k* grows | bootstrap | 10% of the slope, at least .05 per log *n* (a direction: relative rule) | measured |

Also shown: Stan's divergent transitions.

**Choices.** *n* from 50 to 5,000 for the rates: a slope against log *n* needs a wide range. *k* at most 90 for growing *k*: a model has *k* + 3 parameters, and bridge sampling has been shown to work up to 96 (Gronau et al., 2017). Not included: whether information consistency starts exactly at *n* = *k* + 3 (it needs *R*² within 10⁻⁸ of 1, beyond what a sampler resolves).

------------------------------------------------------------------------

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria for Bayesian model choice with application to variable selection. *Annals of Statistics*, *40*(3), 1550–1577.
- Burton, A., Altman, D. G., Royston, P., & Holder, R. L. (2006). The design of simulation studies in medical statistics. *Statistics in Medicine*, *25*(24), 4279–4292.
- Chow, Y. S., & Robbins, H. (1965). On the asymptotic theory of fixed-width sequential confidence intervals for the mean. *Annals of Mathematical Statistics*, *36*(2), 457–462.
- Flegal, J. M., & Gong, L. (2015). Relative fixed-width stopping rules for Markov chain Monte Carlo simulations. *Statistica Sinica*, *25*(2), 655–675.
- Glynn, P. W., & Whitt, W. (1992). The asymptotic validity of sequential stopping rules for stochastic simulations. *Annals of Applied Probability*, *2*(1), 180–198.
- Gronau, Q. F., et al. (2017). A tutorial on bridge sampling. *Journal of Mathematical Psychology*, *81*, 80–97.
- Gronau, Q. F., Singmann, H., & Wagenmakers, E.-J. (2020). bridgesampling: An R package for estimating normalizing constants. *Journal of Statistical Software*, *92*(10), 1–29.
- Jones, G. L., Haran, M., Caffo, B. S., & Neath, R. (2006). Fixed-width output analysis for Markov chain Monte Carlo. *Journal of the American Statistical Association*, *101*(476), 1537–1547.
- Koehler, E., Brown, E., & Haneuse, S. J.-P. A. (2009). On the assessment of Monte Carlo error in simulation-based statistical analyses. *The American Statistician*, *63*(2), 155–162.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11), 2074–2102.
- Pawel, S., Kook, L., & Reeve, K. (2024). Pitfalls and potentials in simulation studies: Questionable research practices in comparative simulation studies allow for spurious claims of superiority of any method. *Biometrical Journal*, *66*(1), e2200091.
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S. (2023). Workflow techniques for the robust use of Bayes factors. *Psychological Methods*, *28*(6), 1404–1426.
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., & Pawel, S. (2024). Simulation studies for methodological research in psychology: A standardized template for planning, preregistration, and reporting. *Psychological Methods*. Advance online publication.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788.
