# flatcor-validation: do flatcor's three priors work?

The validation of [flatcor](https://github.com/dm-k1/flatcor), an R package
for objective Bayes factors and estimation in linear regression: whether each
of its three default priors works correctly within the Bayesian framework.
The argument for the priors, the priors and their derivations, and the paper
are in [Objective-bayes](https://github.com/dm-k1/Objective-bayes).

## The three priors

| | What is flat | Implied $R^2$ | One coefficient's prior |
|---|---|---|---|
| **flat partial r** (`prior = "flat_partial_r"`) | every valid correlation matrix of $(y, x_1 \ldots x_k)$ | Beta($k/2$, 1), grows with $k$ | a flat partial correlation, for every $k$ |
| **flat R** (`prior = "flat_r"`) | the multiple correlation $R$ | Beta(½, 1), for every $k$ | its share of the flat total; narrows as $k$ grows |
| **flat R²** (`prior = "flat_r_squared"`) | $R^2$ | Beta(1, 1), for every $k$ | its share of the flat total; narrows as $k$ grows |

All three are $R^2$ ~ Beta($a$, 1) with $a$ = $k/2$, ½ or 1, and divide the
total over the predictors without preference. As a density on the
coefficients:

```
P(beta) = P(total) x P(this combination | total) x unit conversion
          R^2 ~ Beta(a, 1)     every split equally likely   noise units -> data units
```

flatcor does everything with one Stan model, the same for the three priors:
it draws the prior, gives the posterior, and gives the Bayes factor between
any two linear models, from each model's marginal likelihood estimated by
bridge sampling. Why these three, and why the default is flat partial r: the
paper's repository.

## The validation

What someone needs to trust that each prior works within the Bayesian
framework: three things for every prior, each checked in two ways, with
checks the field already uses.

| | 1. Prior sampling | 2. Posterior estimation | 3. Bayes factors |
|---|---|---|---|
| **a. Computed correctly?** | 1a: the prior draws against the exact prior | 2a: simulation-based calibration | 3a: calibration of the Bayes factor; bridge sampling's error |
| **b. Behaves like a default?** | 1b: invariance | 2b: estimation against a known truth | 3b: the Bayarri et al. (2012) criteria |

Every study reports its estimates against what is expected, with their Monte
Carlo errors; nothing is judged pass or fail per row.

> **TO REFRESH from the full run:** a short summary of the six studies'
> results for each prior, with the figure that shows them best.

## Reading order

Background, in [Objective-bayes](https://github.com/dm-k1/Objective-bayes):
[the three priors](https://github.com/dm-k1/Objective-bayes/blob/main/docs/2_the_three_priors.md)
(part 2: the family $R^2 \sim$ Beta($a$, 1), what each keeps flat and
expects, how the total is divided, what the package computes) and
[their derivations](https://github.com/dm-k1/Objective-bayes/blob/main/docs/appendix/A_derivations.md)
(appendix A: the density factor by factor, the Stan model, the Bayes factors
and the g-prior forms behind them, the default-prior criteria). The table
above is enough to follow the validation. In this repository:

| | Document | What it is |
|---|---|---|
| 1 | [`docs/validation_tests.md`](docs/validation_tests.md) | **What the validation tests**: per study its design, what is tested, what is expected and why, with what Monte Carlo error and how many datasets, and the data sets |
| 2 | [`validation.ipynb`](validation.ipynb) | **The results**: per study its question, its estimates against what is expected, and figures |

Working notes: [`docs/working/`](docs/working/) (the plan of the redesign, the
run-time estimate).

## Names used

| Name used | Means | Also appears as |
|---|---|---|
| **total** | the explained variance $R^2$ of the whole model (flat R keeps its square root $R$ flat) | "shared total" |
| **shape** $a$ | the first parameter of $R^2 \sim$ Beta($a$, 1): $k/2$, ½ or 1 | `r_squared_shape` (Stan), `prior_shape()` |
| **share** | the fraction of the total that one predictor gets; Beta(½, (k−1)/2) for every $a$ | $c^2$ (squared alignment, appendix A section 7), $B$ (appendix A section 8) |
| **noise units**, $f$ | an effect divided by the noise SD $\sigma$; $f^2$ is Cohen's $f^2$ | |
| **Bayes factor for a term** | the model with the term against the model without it, each with its own prior | |
| **references** | JZS default (Cauchy, scale 0.354) and Zellner–Siow (scale 1) | |

## Using the priors: what to watch out for

- **Choosing the prior.** Without `prior =`, flatcor uses flat partial r and
  says so once per session. What each prior assumes:
  [`2_the_three_priors.md`](https://github.com/dm-k1/Objective-bayes/blob/main/docs/2_the_three_priors.md) (in Objective-bayes).
- **Shrinkage at small $n$.** With several predictors and a small sample, the
  prior pulls estimates towards zero (study 2b shows how much).
- **Leaving the null model.** Against the intercept-only model, evidence for
  the null grows like $a \log n$: ½ log n for flat R and log n for flat
  R-squared, whatever the number of predictors, and $\tfrac k2 \log n$ for
  flat partial r.
- **Under flat partial r the claimed $R^2$ grows with $k$.** The median is .63
  with 3 predictors and .87 with 10, so a large model with a modest fit is
  penalised.
- **Flat R-squared with one predictor is not Jeffreys' prior.** Its prior on
  $r$ has density $|r|$, which gives little weight to correlations near zero.
- **Information consistency needs $n \ge k + 3$**, for all three priors.
- **Bayes factors are estimates.** Each comes with a Monte Carlo error
  (bridge sampling's approximate error, which can be optimistic), and models
  up to about 90 predictors stay within the sizes bridge sampling has been
  shown to work for.
- **Scope.** Linear regression with a normal outcome, an intercept and fixed
  effects (continuous predictors, factors tested as one block, interactions,
  polynomials, splines). Not supported: group-level effects, other families,
  distributional models, models without an intercept.

## Running the validation

Restore the environment with `renv::restore()` (it installs flatcor at the
pinned version, with targets, crew and the other packages), or use the
container image of `cluster/Dockerfile`, published by GitHub Actions as
`ghcr.io/dm-k1/flatcor-validation`. Then:

```bash
VALIDATION_SMOKE=1 Rscript run_all.R   # first: a corner of every study, about 10 minutes
Rscript run_all.R                      # the studies, then the figures
VALIDATION_FRESH=1 Rscript run_all.R   # the same, from scratch
```

The smoke test runs every step and every figure on a small part of each
design, from scratch, with its results in `cache/` (never mixed with real
ones).

or open [`validation.ipynb`](validation.ipynb), which runs through
`run_all.R` too. A stopped run continues where it was. It refuses to
continue with another flatcor than the one it started with (each results
file records the flatcor commit), so one results file never mixes two
versions; start from scratch then. (`targets::tar_make()` runs the pipeline
without that check: for working on a single step.)
[`cluster/README.md`](cluster/README.md) describes the container image, how
long the full run takes, and how to run it on a cluster (Apptainer, a SLURM
job, resubmitting until it is done).

Every push to `main` is checked automatically
([`.github/workflows/check.yml`](.github/workflows/check.yml)): lintr, the
links, the workflow files, and the smoke test inside the published image.

The smoke test is meant to run the way a real run does, because a bug can
hide anywhere the two differ. So it runs in the same image (the same
packages and CmdStan) and with several workers (3, by pretending to have 16
cores), since workers clashing over shared files only happens with more
than one. Two differences remain, on purpose: it stops after two rounds
(round three and later run the same code as round two), and it runs a
corner of each design.

## What's where

| Path | What it is |
|---|---|
| `_targets.R` | The pipeline: the steps of every study, in order |
| `run_all.R` | Runs the pipeline, then the figures |
| `R/shared/` | Used by every step: the random streams, the pipeline's helpers, the reporting rule, the results file, the comparison methods |
| `R/1_data/` | One file per data set, numbered in the order the studies first use them (what each is and which study needs it: `docs/validation_tests.md`, "The data sets") |
| `R/2_studies/` | The six studies, one folder per cell, each with `design.R` and `code.R` |
| `R/3_figures/` | One script per study that draws its figures from its results file |
| `R/4_notebook/` | The notebook's helpers: they present the results, they compute nothing |
| `R/timing/` | The run-time estimate: time a study on a fraction of its datasets, scale up to the full run |
| `results/` | Each study's tables and setup, `studies/study_<cell>.rds`, and its figures, `figures/<cell>_*.png` |
| `validation.ipynb` | The notebook (R kernel) |
| `docs/` | What the validation tests, working notes |
| `cluster/` | The container image, the run-time estimate, and how to run the validation on a cluster |
| `.github/workflows/` | The checks on every push (`check.yml`), the image build (`image.yml`) and the run-time estimate, on demand (`timing.yml`) |
| `renv.lock`, `.lintr` | The pinned packages; the code style |

## Author

Diemithry M. Kloppenburg, University of Amsterdam.

## License

MIT (`LICENSE.md`).
