# Validation redesign

Working plan for the branch `validation-redesign`: what the validation is,
the decisions behind it, and what is left. When the redesign is done, this
file is replaced by the finished validation document. The plans before 14
November 2026 (the ten-study structure, the run of 4 November 2026 and its
documents in `docs/working/old_run/`) are in the git history, last in commit
9ffd1bf.

## What the validation is for

What someone needs to trust that each of the three priors works within the
Bayesian framework, and nothing else:

| Study | Question | Standard practice | Data |
|---|---|---|---|
| 1 | Does the code draw the prior it claims? | prior check against the exact distribution | none |
| 2 | Are the Bayes factors computed correctly? | against routes that share no code with the package: numerical integration of a second derivation, a 200-bit ground truth, brute force, bridge sampling on the Stan model | none / few |
| 3 | Are the posterior and the Bayes factor calibrated? | SBC with its divergence diagnostics (Talts et al., 2018); Bayes factor calibration (Schad et al., 2023), on data from each prior | B |
| 4 | Does each prior meet the objective-Bayes criteria? | Bayarri, Berger, Forte & García-Donato (2012): predictive matching, information consistency, intrinsic consistency, model-selection consistency (fixed and growing k), invariance | none / D |

How the priors behave in use (how often the test of a term finds an effect,
confirms a null or stays undecided as k and n change; separating a null from
a small effect; estimation against OLS) is a separate study, for later. Its
code, with data set A and its factor extension, is in commit 9ffd1bf
(`R/studies/05_behaviour.R`).

## Decisions

1. **Three ideas, each in the standard Bayesian sense:** the Bayes factor
   between two models (each with its own prior); a term as the model with it
   against the model without it (a factor is the same); the posterior from
   flatcor's Stan model.
2. **One implementation for all three priors.** Each prior is R² ~ Beta(a, 1)
   with a uniform direction, and only a changes (flat partial r k/2, flat R
   1/2, flat R² 1). There are two implementations in all: the Bayes factor of
   a model against the intercept-only model (the likelihood ratio averaged
   over the prior by `integrate()`) and the Stan model. Every check runs the
   same way for every prior; a route that works for some priors and not
   others is used for none.
3. **The conditional test and `bayes_factor_per_term()` are gone**, from the
   validation and from flatcor (2.0.0). The support curves stay: they use no
   prior.
4. **Only what is needed to trust the prior** (above). Study 4 stays because
   "objective default" is part of what flatcor offers. Studies on coherence
   (Savage–Dickey), variable selection and behaviour are out.
5. **Two reference methods, JZS and Zellner–Siow** (BayesFactor), in study 4
   next to the three priors, and later in the behaviour study.
6. **Established frameworks, all in R:** `SimDesign` for the simulations,
   `SBC`, `bridgesampling`, `Rmpfr`. (`BAS` and `rsimsum` left `renv.lock`
   with the studies that used them; the behaviour study can add `rsimsum`
   back.)
7. **The number of datasets follows from a target Monte Carlo SE, fixed in
   advance** (Morris et al., 2019); every reported number carries its Monte
   Carlo SE, and every check's threshold is relative to it.
8. **Datasets are defined by their ID**: each has its own random stream, so
   results do not depend on the machine, the order or the number of cores.
9. **The package lives in its own repository**,
   [dm-k1/flatcor](https://github.com/dm-k1/flatcor), and this repository
   validates a fixed release of it: v2.0.0 (commit cec149c), pinned in
   `renv.lock`. Internal functions some checks need are taken from its
   namespace in `R/shared/settings.R`.

## Status (14 November 2026)

- flatcor 2.0.0 released and pinned. One route for every prior, checked
  against `integrate()` of the g-prior form (at most 1.5e-7) and a nested
  third route; the old closed forms were the less accurate side.
- The four studies are written and run end to end in smoke mode, through
  `run_all.R` and the notebook; appendix C describes them.
- The container image for flatcor 2.0.0 and a new run-time estimate run on
  GitHub Actions (workflows "Validation image" and "Run-time estimate").
  **TO REFRESH:** `docs/working/run_time_estimate.md` once the estimate for
  the four studies is in.

## Still to do

1. The run-time estimate for the four studies (above).
2. The full run of studies 1 to 4, each where it fits (`cluster/README.md`).
3. Rewrite the docs from the results: `docs/3_results.md` (new), the README's
   findings, and the placeholders marked **TO REFRESH** in the README,
   `docs/2_the_three_priors.md` and appendices A and B.
4. Not yet tested: automated tests for seeding and resuming (seeding is
   confirmed by hand: the same check lines on 1 and 4 cores); the SLURM array
   and the Apptainer commands on a real cluster (only mock SLURM variables so
   far).

For a flatcor 2.0.1, found while porting, not yet released:

- With zero residual degrees of freedom the likelihood ratio drops the point
  v = 0 (0 × log 0), so check 4.1 gives 1.0033 instead of 1 at k = 1. The
  package's own Bayes factors never reach this (they return 1 there).
- The compiled Stan model is kept in one folder for every version, so after
  a reinstall a stale executable can crash (seen once; a recompile fixed it).
- Stan's "proposal is about to be rejected" messages in a chain's first few
  iterations (a step that takes R to 1 or sigma to 0; no divergent
  transitions) could be silenced.

## Pitfalls already handled

SimDesign masks a function named `summarise` and reserves all-capital names;
every `runSimulation()` has its own resume file (`save_details`);
`SimExtract()` needs a `summarise` function; Stan models are compiled into
`cache/` (or a temporary folder), never into the installed package, which is
read-only in an Apptainer image; each study caches per flatcor build
(`cache/study_NN/flatcor_<build>/`), so results from another build are never
reused.

## Layout

```
R/
  shared/           settings, precision targets, comparison methods
  data/             the shared data sets B and D, and the ID-based streams
  studies/          the four studies, 01-04
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
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11),
  2074–2102. https://doi.org/10.1002/sim.8086
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S.
  (2023). Workflow techniques for the robust use of Bayes factors.
  *Psychological Methods*, *28*(6), 1404–1426.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
  Validating Bayesian inference algorithms with simulation-based
  calibration. arXiv:1804.06788.
