# Running the validation elsewhere, and how long it takes

The validation is one pipeline ([`_targets.R`](../_targets.R)); `Rscript
run_all.R` runs it and then the figures. It saves every finished piece as it
goes, so it can run anywhere a run may stop (a laptop, Colab, a cloud
machine) and continue where it was, as long as the project folder (with its
`_targets/` store) survives. This folder holds the two things that help with
that: a container image with the exact environment, and the run-time
estimate.

## The container image

`cluster/Dockerfile` builds R 4.6.1, the packages of `renv.lock` and CmdStan
2.40.0. GitHub Actions ([`.github/workflows/image.yml`](../.github/workflows/image.yml))
builds it whenever the environment changes, runs the whole pipeline in smoke
mode inside it, and only then publishes it as
`ghcr.io/dm-k1/flatcor-validation:<branch>` (`latest` on `main`). To run the
validation in it, from the project folder:

```bash
docker run --rm -v "$PWD":/project ghcr.io/dm-k1/flatcor-validation:study-rework
```

Without the image, `renv::restore()` installs the same packages; CmdStan comes
from `cmdstanr::install_cmdstan()`.

## How long it takes

The workflow "Run-time estimate" ([`.github/workflows/timing.yml`](../.github/workflows/timing.yml))
measures it on GitHub's runners, for free on a public repository. It starts
when the studies' code changes on `validation-redesign` or `study-rework`
(a newer push cancels a run still going), and from Actions > Run-time estimate >
Run workflow once it is on `main`. Each study runs its full design on 2% and
5% of its datasets, each on its own 4-core runner; `cluster/estimate_run_time.R`
separates the fixed part of a study's cost from the part per dataset and reads
off the full run. The run's summary page shows the hours per study on 4, 8, 16
and 32 cores.

The same can be run anywhere, a study at a time:

```bash
VALIDATION_TIMING=0.02 Rscript cluster/time_study.R 3a
VALIDATION_TIMING=0.05 Rscript cluster/time_study.R 3a
Rscript cluster/estimate_run_time.R
```

Timing runs write only to `cache/timing/`, so they never mix with a real run.

A SLURM array script for university clusters was here until 17 November 2026
(in the git history); the pipeline can run on such a cluster through crew's
cluster controllers if that is ever needed.
