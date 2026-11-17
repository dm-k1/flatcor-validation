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

Each cell is one study folder in `studies/`, named after its heading (1-3)
and row (a: computed correctly; b: behaves like a default), so the folder
reads like the table:

```
studies/
  1a_prior_sampling/      prior, computed correctly: Stan's prior-only draws = the exact prior
  1b_invariance/          prior, behaves like a default: rescaling or recoding changes nothing
  2a_posterior/           posterior, computed correctly: SBC and Stan's diagnostics
  2b_estimation/          posterior, behaves like a default: bias, coverage, prediction against OLS
  3a_bayes_factors/       Bayes factor, computed correctly: calibration, bridge sampling's error
  3b_default_criteria/    Bayes factor, behaves like a default: the Bayarri criteria
```

Each folder holds `design.R` (the design: cells, numbers of datasets) and
`code.R` (the steps: generate, analyse, summarise, save). What each study
tests, what it expects and why, with what Monte Carlo error, how its data are
made and the reasons for every choice are in one document,
[`docs/validation_tests.md`](../validation_tests.md). The pipeline
[`_targets.R`](../../_targets.R) runs the studies (targets): every piece of
work is saved as it finishes, a stopped run continues where it was, and only
what is not done yet or has changed reruns. Each study writes its tables to
`results/studies/study_<cell>.rds`; the figures are a separate step, one
script per study in `figures/`, that reads those files. `Rscript run_all.R`
runs the pipeline and then the figures. Items are numbered within their
study (3b.1, 3b.2, ...), so a number says which cell it belongs to.
Invariance sits under the prior because it holds by how the prior is built,
although it is tested on the Bayes factor. The notebook follows the same
structure.

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
- Every study reports this way. No study has a tripwire yet; the notebook
  shows any that a study records.

How the priors behave in use beyond estimation (how often the test of a
term finds an effect, confirms a null or stays undecided) is a separate
study, for later. Its earlier code, with a factor extension of the data, is
in commit 9ffd1bf (`R/studies/05_behaviour.R` there, parts A and B).

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

### What the bridgesampling authors say (read 17 November 2026)

From Gronau et al. (2017), "A tutorial on bridge sampling", and Gronau,
Singmann & Wagenmakers (2020), the package paper (its Discussion lists four
points of practical advice). What flatcor and the studies take from it:

| Point | What the authors say | What we do |
|---|---|---|
| Number of parameters | No stated limit. Bridge sampling "performs accurately even in high-dimensional parameter spaces" such as hierarchical models (2020, introduction). The normal proposal "might produce unstable estimates" for high-dimensional posteriors that are clearly not multivariate normal (2017). Their largest worked example has 96 parameters (2017: a hierarchical model, 30 participants x 3 parameters + 6 group-level), with 60,000 posterior draws; the package paper's examples are smaller (a GLMM with 31 random effects; factor models) | Our model has k + 3 parameters, so k up to about 90 stays within what the authors have shown. Beyond that we make no claim |
| Draws needed | Testing needs "about an order of magnitude more posterior samples than estimation" (2020, rule of thumb in the first example and the Discussion); how many depends on the model and the accuracy wanted | 10 times the draws used for estimation: 4 chains of 10,000 instead of 1,000 |
| Normal or Warp-III | Warp-III is always at least as precise as the normal method and takes about twice as long; it helps when the joint posterior (on the unconstrained scale) is skewed, and they suggest trying it "when in doubt" (2020, sections 2.1-2.2 and the Discussion). With several modes neither is sure to work | The normal method by default; Warp-III where the posterior is skewed (a candidate cause of the flat R-squared misses at k = 1, open question below) |
| The error estimate | `error_measures()` gives an approximate relative mean-squared error (Frühwirth-Schnatter, 2004) that allows for autocorrelated MCMC draws; available for the normal method only, as for Warp-III it was unreliable in their experience. `repetitions` > 1 shows the variability from the proposal draws only, with the posterior draws fixed. Their "gold standard": repeat everything with fresh posterior draws and look at the spread of the Bayes factor itself (2020, section 2.3 and the Discussion) | The approximate error for every Bayes factor; repeated runs with fresh posterior draws in study 3a, to check that error |
| Requirements | All normalising constants kept (in Stan `target +=` rather than `~`); no discrete parameters (marginalise them out); converged chains with enough representative draws; parameters on the real line (the package transforms bounded ones) | flatcor's model uses `target +=` throughout; the sign of R with one predictor is marginalised (`log_mix`); Stan's diagnostics are reported |
| Improper priors | Not addressed in either paper | Our condition above (the same improper priors on intercept and noise in both models) is ours to state, with the standard argument (Jeffreys; Bayarri et al., 2012) |
| The prior | Bayes factors depend on it; they point to default priors with properties such as model-selection and information consistency, citing Bayarri et al. (2012) | Study 3b checks exactly those properties |

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
   and the validation says where that is. For bridge sampling that is up to
   about 96 parameters, the largest case its authors have shown (above), so
   k up to about 90; growing k (Bayarri's model-selection consistency as k
   grows with n) is tested within that range, not at k = 500.
6. **Two reference methods, JZS and Zellner–Siow** (BayesFactor), next to
   the three priors in the default-prior checks and in study 2b (with OLS
   there), and later in the behaviour study.
7. **Established frameworks, all in R:** `targets` (with `crew` for parallel
   workers) runs the studies, `SBC`, `bridgesampling`, `cmdstanr`. Not
   SimDesign: it loads RcppParallel, whose TBB library breaks cmdstanr's
   log density, which bridge sampling needs; and targets keeps the steps
   readable, saves every piece as it finishes and resumes a stopped run.
   One backend for all Stan work: cmdstanr (rstan is about 1.6 times slower
   per Bayes factor here; flatcor supports both through its `backend`
   argument).
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

- **flatcor 3.0** (in development) is on flatcor's `main` (commit 9a80b87)
  and pinned in `renv.lock`: one Stan model for the three priors, with the
  intercept-only model (K = 0) and every normalising constant;
  `bayes_factor_models()` by bridgesampling (4 chains of 10,000 draws per
  model); the formula Bayes factors, the likelihood ratio code and the
  support curves removed; the compiled model in a folder per version. Its
  tests check the intercept-only model against its closed form and each
  prior against an independent reference (g-mixture or direction form).
  The Kummer-series route stays in a git stash on flatcor, not released.
- cmdstanr can only give a model's log density (which bridgesampling needs)
  in the R session that compiled it, so flatcor compiles that version once
  per session, into a temporary folder (commit 9a80b87, on `main`).
- **The studies follow the six cells** as a targets pipeline (branch
  `targets-3a`, which builds on `study-rework`): `studies/<cell>/` with
  `design.R` and `code.R`, `_targets.R`, figure scripts in `figures/`,
  `run_all.R` as the run-all, and `docs/validation_tests.md` as the one
  document of what is tested. The old scripts in `R/studies/`, the shared
  `settings.R` and `precision.R`, appendix C (now `validation_tests.md`),
  SimDesign and the SLURM array are gone (in the git history); the
  motivating simulation keeps its helpers in `R/docs/helpers.R`.
- **Study 3a shares its data across the priors:** the same predictors, noise
  and model assignment for the three priors, only the effects drawn from
  each prior (common random numbers); intercept 0 and noise SD 1.
- **From the old study 4, dropped:** whether information consistency starts
  exactly at n = k + 3 (it needs R-squared within 10^-8 of 1 at n = k + 2,
  beyond what a sampler resolves; appendix A, section 11, keeps the
  derivation). Growing k now stops at k = 90.
- **Smoke test of the whole pipeline:** **TO REFRESH**.

## Still to do

1. The run-time estimate for the six studies; with it, set the provisional
   numbers of datasets of 3a and 3b (their `design.R`).
2. The full run (`Rscript run_all.R`; `cluster/README.md`).
3. The docs from the results: `docs/3_results.md` (new), the README's
   findings, and the places marked **TO REFRESH** in the README,
   `docs/2_the_three_priors.md`, appendices A and B,
   `docs/validation_tests.md` and `docs/working/run_time_estimate.md`.

**Open questions:**

- Resolved (17 November 2026): the two bridge-sampling misses for flat R²
  with one predictor (timing run of 248d5d6; an earlier smoke run) were not
  bias. 50 runs of that case (n = 80, 16,000 draws each) averaged -1.7303
  (MCSE 0.0004) against the exact -1.7301, but spread with SD 0.0025 where
  bridgesampling reported 0.0018: its approximate error is about 1.4 times
  too small there, so 3 of 50 runs lay beyond 3 reported errors. The
  repeated runs with fresh posterior draws in 3a (the authors' gold
  standard) give the error that is reported; without per-row pass/fail this
  cannot become a false failure. Only this case was run.
- The cost of Bayes factor calibration (3a) and the Bayarri checks (3b) with
  bridge sampling: each Bayes factor is now two Stan fits with 10 times the
  estimation draws, where it was one formula (about 40 s for n = 60, k = 3
  on this Mac). The numbers are now provisional (in each study's `design.R`):
  1,000 datasets per cell for calibration, 100 per setting and n for the
  rates, 50 per setting for growing k. The run-time estimate sets them.

## Pitfalls already handled

cmdstanr gives a model's log density (for bridgesampling) only in the R
session that compiled it, so flatcor compiles that version per session into
a temporary folder; that compile fails once RcppParallel is loaded (SimDesign
and rstan load it), so the pipeline loads neither. Stan models are compiled
into `cache/` or a temporary folder, never into the installed package, which
is read-only in a container. The pipeline's store (`_targets/`) records which
version of the code and of flatcor each result came from, and reruns what is
outdated.

## Code style

Every R file follows the tidyverse style (styler) and passes lintr with
the project's `.lintr` (lines up to 100 characters; names may end in a
capital symbol, as `log_BF10`; "no visible binding" is off,
because the pipeline sources the files together). Within a function, a
blank line and a short comment separate each step; longer calls put one
argument per line. Check with `lintr::lint_dir()`. Run styler in a UTF-8
locale (`LC_ALL=en_US.UTF-8`): in an ASCII locale it writes characters such
as ² as `<U+00B2>`.

## Layout

```
_targets.R          the pipeline: the steps of every study, in order
run_all.R           runs the pipeline, then the figures
studies/            one folder per study (1a-3b): design.R and code.R
R/
  shared/           pipeline helpers, the reporting rule, the results file, comparison methods
  data/             one file per data set (prior-only designs, example data, data from the
                    prior, known slopes, set R², nested models, growing k) and the
                    building blocks (ID-based random streams); docs/validation_tests.md
  notebook/         present.R, the notebook's helpers (they compute nothing)
  docs/             the motivating simulation, its helpers, and the figures of the docs
figures/            one script per study, and the figures they draw (tracked)
results/            studies/study_<cell>.rds per study: its tables and the setup of the run (tracked)
validation.ipynb    the notebook (R kernel): per study its question, tables and figures
cluster/            the container image (Dockerfile) and the run-time estimate
docs/               parts 1 and 2, validation_tests.md, appendix/ (A, B), working/ (this plan)
_targets/, cache/   the pipeline's store, smoke and timing runs (gitignored)
```

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012).
  Criteria for Bayesian model choice with application to variable
  selection. *The Annals of Statistics*, *40*(3), 1550–1577.
- Landau, W. M. (2021). The targets R package: A dynamic Make-like
  function-oriented pipeline toolkit for reproducibility and
  high-performance computing. *Journal of Open Source Software*, *6*(57),
  2959.
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
