# Validation redesign

Working plan for the branch `validation-redesign`: what the validation is,
the decisions behind it, and what is left. When the redesign is done, this
file is replaced by the finished validation document. Earlier plans are in
the git history: the ten-study structure and the run of 4 November 2026 last
in commit 9ffd1bf, the four-study plan built on flatcor's own Bayes factor
formulas in commit 248d5d6.

## What the validation is for

What someone needs to trust that each of the three priors works within the
Bayesian framework, and nothing else. Three things are validated for every
prior, each with checks the field already uses:

1. **Prior sampling**: the software draws the prior it claims.
2. **Posterior estimation**: the posterior is computed correctly.
3. **Bayes factors**: the Bayes factor between two models is computed
   correctly.

For each of the three, two questions: is it computed correctly, and does it
behave the way a default prior should?

| | 1. Prior sampling | 2. Posterior estimation | 3. Bayes factors |
|---|---|---|---|
| **Computed correctly?** | Stan's prior-only draws against the exact prior: R ~ Beta(2a, 1), direction uniform (prior check) | SBC (Talts et al., 2018; `SBC` package, Modrák et al., 2023), with Stan's diagnostics: divergent transitions, R-hat, effective sample size (Vehtari et al., 2021) | bridge sampling with its own error estimate and repeated runs (Gronau et al., 2017, 2020); Bayes factor calibration on data from each prior (Schad et al., 2023) |
| **Behaves like a default?** | invariance: rescaling or recoding the predictors changes nothing (Bayarri et al., 2012) | behaviour in estimation: bias and coverage of a slope, prediction error against OLS, on data with a known truth | predictive matching, information consistency, intrinsic consistency, model-selection consistency (Bayarri et al., 2012) |

Each cell is one study file in `R/studies/`, named after its heading (1-3)
and row (a: computed correctly; b: behaves like a default), so the folder
reads like the table and each file answers one question:

```
R/studies/
  1a_prior_sampling.R      prior, computed correctly: Stan's prior-only draws = the exact prior
  1b_invariance.R          prior, behaves like a default: rescaling or recoding changes nothing
  2a_posterior.R           posterior, computed correctly: SBC and Stan's diagnostics
  2b_estimation.R          posterior, behaves like a default: bias, coverage, prediction against OLS
  3a_bayes_factors.R       Bayes factor, computed correctly: bridge sampling's error, calibration
  3b_default_criteria.R    Bayes factor, behaves like a default: the Bayarri criteria
```

Checks are numbered within their file (3b.1, 3b.2, ...), so a check's number
says which cell it belongs to. The cells differ in cost (1a and 1b take
minutes, 2a, 3a and 3b hours), and as separate files they run as separate
cluster jobs. Invariance sits under the prior because it holds by how the
prior is built, although it is tested on the posterior and the Bayes
factor. The notebook follows the same structure: three sections, each with
a "computed correctly" and a "behaves like a default" part.

No diagnostics of our own: every check is one of the standards above. What
is ours is the choice of settings (n, k, R²), the expected values and the
number of datasets, which follows from a target Monte Carlo SE fixed in
advance (Morris et al., 2019).

**Reporting: estimates against expectations, not pass/fail.** A simulation
study is estimation (Morris et al., 2019), so every cell reports each
quantity as one row:

| Quantity | Expected (fixed in advance) | Estimate | Monte Carlo SE | Distance in MCSEs |
|---|---|---|---|---|

with, under each table, how many estimates lie beyond 2 and 3 Monte Carlo
SEs against what chance gives (about 5% and 0.3%). The reader judges
whether the estimates are sensible; the software gives no verdict per row.

- What is committed in advance is the expected value and the Monte Carlo SE
  target, so a result cannot be judged "fine" after the fact.
- The kind of expectation differs by row. Computed correctly (a): an exact
  value from theory (uniform SBC ranks, a calibration slope of 1, the
  prior's quantiles). Bayarri criteria (3b): often a direction or a limit
  ("grows with n", "equals 1"), reported next to the quantity, such as a
  slope with its Monte Carlo SE against the slope theory gives. Estimation
  (2b): references that are not guaranteed (coverage .95, prediction ratio
  to OLS 1), so the report is descriptive.
- Pass/fail stays only as a tripwire for broken code: exact identities
  (predictive matching gives exactly 1, invariance gives exactly equal
  results), a study that stops with an error, and Stan's diagnostics are
  reported as counts. Per-row pass/fail on Monte Carlo estimates misleads:
  with dozens of estimates some lie beyond 2 MCSE by chance (as check 4.8
  at few draws, and a bridge-sampling difference of 0.006 against a
  threshold of 0.0057, did in the timing run of 248d5d6).
- Study 2b already reports this way; the other cells change with the
  rework, and `record_check()`/`report_checks()` and the notebook's checks
  table shrink to the tripwires.

How the priors behave in use beyond estimation (how often the test of a
term finds an effect, confirms a null or stays undecided) is a separate
study, for later. Its earlier code, with data set A's factor extension, is
in commit 9ffd1bf (`R/studies/05_behaviour.R`, parts A and B).

## The software: flatcor

flatcor does the three things for all three priors in one integration:

| | How | Whose code |
|---|---|---|
| The prior | flatcor's Stan model: R ~ Beta(2a, 1) on R (not R²), direction uniform; only a changes (flat partial r k/2, flat R 1/2, flat R² 1) | ours: stating the prior |
| Prior sampling | the same model with the data switched off (`prior_only = 1`) | Stan |
| Posterior | the same model with the data | Stan |
| Bayes factor | bridgesampling on two fits of the model: with the term and without it | bridgesampling, with a short wrapper of ours |

Sampling R instead of R² keeps every prior's density finite at zero (flat R
is flat on R; on R² it would be infinite at 0), which is why SBC passes
without divergent transitions.

The intercept and the noise get improper priors (flat and 1/sigma). The
Bayes factor is still valid because both models give them exactly the same
prior, so the arbitrary constant cancels; the documentation states this as a
condition.

## Decisions

1. **Three things, each in the standard Bayesian sense:** the prior, the
   posterior from Stan, and the Bayes factor between two models (a term
   tested as the model with it against the model without it; a factor is the
   same).
2. **One implementation for all three priors:** one Stan model in which only
   a changes. Every check runs the same way for every prior.
3. **No numerical software of our own.** Exact Bayes factor formulas
   (flatcor 2.0.0's likelihood ratio averaged by `integrate()`, and the
   Kummer-series route tried on 17 November 2026) kept breaking in corners:
   many predictors (k of about 90 or more), near-perfect fits with tiny n,
   R² = .9999 with n = 10⁵. Fixing those is numerical work, not validation of
   a prior. The Bayes factor comes from bridgesampling, an established
   package. A known limit of Stan or bridgesampling is stated, not solved.
4. **Only what is needed to trust the prior** (above). Coherence
   (Savage–Dickey), variable selection and behaviour are out.
5. **The range is stated.** Checks run where the tools are known to work,
   and the validation says where that is. Growing k (Bayarri's
   model-selection consistency as k grows with n) is tested at moderate k,
   not at k = 500, because bridge sampling needs far more draws as the
   number of parameters grows.
6. **Two reference methods, JZS and Zellner–Siow** (BayesFactor), next to
   the three priors in the default-prior checks and in study 2b (with OLS
   there), and later in the behaviour study.
7. **Established frameworks, all in R:** `SimDesign` for the simulations,
   `SBC`, `bridgesampling`, `cmdstanr`.
8. **The number of datasets follows from a target Monte Carlo SE, fixed in
   advance** (Morris et al., 2019); every reported number carries its Monte
   Carlo SE and is shown against its expected value (the reporting rule
   above), not judged pass or fail.
9. **Datasets are defined by their ID**: each has its own random stream, so
   results do not depend on the machine, the order or the number of cores.
10. **The package lives in its own repository**,
    [dm-k1/flatcor](https://github.com/dm-k1/flatcor), and this repository
    validates a fixed release of it, pinned in `renv.lock`.

## Status (17 November 2026)

- flatcor 2.0.0 is released and pinned. Its Stan model already samples all
  three priors in R; SBC passes for all three (k up to 30) without divergent
  transitions.
- Its Bayes factors are its own formulas, which this plan replaces. The
  timing run of commit 248d5d6 found studies 1-3 passing and study 4 stuck at
  growing k (flat R with about 90 or more predictors: a Bessel function
  underflows).
- The Kummer-series route is on the local flatcor branch `kummer-series`,
  uncommitted. It passes flatcor's 124 tests but fails 15 of study 2's 792
  settings (R² = .9999, n = 10⁵) and runs out of memory on check 4.3. Not to
  be released.
- The studies in `R/studies/` (01-04) and the notebook still follow the
  previous plan; they are reworked below.
- Study 2b is written (`R/studies/2b_estimation.R`, ported from part C of
  the old behaviour study, with data set A in `R/data/data_set_A.R`), passes
  its smoke test, and is in `run_all.R`, the SLURM array (task 5), the timing
  workflow and the notebook. It uses only Stan fits, so it does not wait for
  the flatcor rework.

## Still to do

**flatcor (a new major version: the formula Bayes factors go):**

1. A Bayes factor function on the Stan model: fit the model with and without
   the term, `bridgesampling::bridge_sampler()` on each, return the Bayes
   factor with its error estimate.
2. The intercept-only model (k = 0) for bridgesampling: the Stan model
   currently needs k >= 1.
3. Remove `log_bayes_factor_vs_intercept_only()`, the likelihood ratio code
   (`R/likelihood.R`) and the tests built on them. The support curves use the
   likelihood ratio; decide whether they stay.
4. The two open items from 2.0.0: a separate folder for the compiled Stan
   model per version (a stale executable crashed once after a reinstall);
   optionally silence Stan's "proposal is about to be rejected" messages in a
   chain's first iterations.

**The validation, by cell** (from the current studies 01-04):

| Now | Becomes |
|---|---|
| 01 prior sampling | `1a_prior_sampling.R`, on Stan's prior-only draws |
| 04's invariance check (4.11) | `1b_invariance.R` |
| 03.1 SBC and 03.2 divergences | `2a_posterior.R` |
| 03.3 Bayes factor calibration | `3a_bayes_factors.R`, with bridge sampling, plus bridge sampling's own error and repeated runs |
| 04 without 4.11 | `3b_default_criteria.R`, with bridge sampling; growing k at moderate k |
| 02 Bayes factors against independent routes | gone: it checked our formulas. An independent route for bridge sampling (such as flat partial r's g-integral) is not needed; Bayes factor calibration shows whether the Bayes factors are right end to end |
| (commit 9ffd1bf, `05_behaviour.R`, part C) | `2b_estimation.R`: done |

`run_all.R`, `cluster/slurm_array.sh`, the timing workflow and the notebook
follow the new file names.

**Open questions:**

- How far the bridgesampling authors take it (number of parameters, draws
  needed): look up Gronau et al. (2017, 2020) before fixing the range.
- Bridge sampling missed twice for flat R² with one predictor (the timing run
  of 248d5d6: off by 0.006 against 3 x its error 0.0019; and an earlier
  smoke run): chance or a pattern.
- The cost of Bayes factor calibration with bridge sampling: two fits per
  dataset instead of one formula.

**Then:** the run-time estimate, the full run, and the docs from the
results (`docs/3_results.md`, the README's findings, the placeholders marked
**TO REFRESH** in the README, `docs/2_the_three_priors.md` and appendices A
and B; appendix A's derivation of the likelihood ratio and appendix C's
description of study 2 change with this plan).

## Pitfalls already handled

SimDesign masks a function named `summarise` and reserves all-capital names;
every `runSimulation()` has its own resume file (`save_details`);
`SimExtract()` needs a `summarise` function; Stan models are compiled into
`cache/` (or a temporary folder), never into the installed package, which is
read-only in an Apptainer image; each study caches per flatcor build
(`cache/study_NN/flatcor_<build>/`), so results from another build are never
reused.

## Layout (current; changes with the rework)

```
R/
  shared/           settings, precision targets, comparison methods
  data/             the shared data sets B and D, and the ID-based streams
  studies/          the four studies, 01-04 (previous plan); one file per
                    cell after the rework (1a-3b, above)
  notebook/         present.R, the notebook's helpers (they compute nothing)
  docs/             the motivating simulation and the figures of the docs
results/            studies/study_NN.rds per study (tables, checks, figures,
                    and where, when and how long it ran); session_info/;
                    run_summary.md from run_all.R (tracked)
figures/            every figure of the studies (tracked)
validation.ipynb    the notebook (R kernel): per study its question, checks,
                    tables and figures; reads results/ or reruns the studies
run_all.R           the same without the notebook, one study after another
cluster/            Dockerfile, slurm_array.sh, the timing scripts, README.md
docs/               parts 1 and 2, appendix/ (A-C), working/ (this plan)
cache/              raw results, logs, smoke and timing runs (gitignored)
```

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012).
  Criteria for Bayesian model choice with application to variable
  selection. *The Annals of Statistics*, *40*(3), 1550–1577.
- Chalmers, R. P., & Adkins, M. C. (2020). Writing effective and reliable
  Monte Carlo simulations with the SimDesign package. *The Quantitative
  Methods for Psychology*, *16*(4), 248–280.
- Gronau, Q. F., Sarafoglou, A., Matzke, D., Ly, A., Boehm, U., Marsman, M.,
  Leslie, D. S., Forster, J. J., Wagenmakers, E.-J., & Steingroever, H.
  (2017). A tutorial on bridge sampling. *Journal of Mathematical
  Psychology*, *81*, 80–97.
- Gronau, Q. F., Singmann, H., & Wagenmakers, E.-J. (2020). bridgesampling:
  An R package for estimating normalizing constants. *Journal of
  Statistical Software*, *92*(10), 1–29.
- Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková,
  K., Gelman, A., & Vehtari, A. (2023). Simulation-based calibration
  checking for Bayesian computation: The choice of test quantities shapes
  sensitivity. *Bayesian Analysis*.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11),
  2074–2102. https://doi.org/10.1002/sim.8086
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S.
  (2023). Workflow techniques for the robust use of Bayes factors.
  *Psychological Methods*, *28*(6), 1404–1426.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
  Validating Bayesian inference algorithms with simulation-based
  calibration. arXiv:1804.06788.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
  (2021). Rank-normalization, folding, and localization: An improved R-hat
  for assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718.
