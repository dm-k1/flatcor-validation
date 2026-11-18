# Running the validation elsewhere, and how long it takes

The validation is one pipeline ([`_targets.R`](../_targets.R)); `Rscript
run_all.R` runs it and then the figures. It saves every finished piece as it
goes, so it can run anywhere a run may stop (a laptop, Colab, a cloud
machine) and continue where it was, as long as the project folder (with its
`_targets/` store) survives and flatcor stays the same: `run_all.R` refuses
to continue with another flatcor commit than the run started with.
`VALIDATION_FRESH=1 Rscript run_all.R` starts from scratch. This folder
holds what helps with that: a container image with the exact environment,
the run-time estimate, and how to run it on a cluster.

## The container image

`cluster/Dockerfile` builds R 4.6.1, the packages of `renv.lock` and CmdStan
2.40.0. GitHub Actions ([`.github/workflows/image.yml`](../.github/workflows/image.yml))
builds it whenever the environment changes, runs the whole pipeline in smoke
mode inside it, and only then publishes it as
`ghcr.io/dm-k1/flatcor-validation:<branch>` (`latest` on `main`). To run the
validation in it, from the project folder:

```bash
docker run --rm -v "$PWD":/project ghcr.io/dm-k1/flatcor-validation:latest
```

Without the image, `renv::restore()` installs the same packages; CmdStan comes
from `cmdstanr::install_cmdstan()`.

## How long it takes

The workflow "Run-time estimate" ([`.github/workflows/timing.yml`](../.github/workflows/timing.yml))
measures it on GitHub's runners, for free on a public repository. It runs on
demand only, from Actions > Run-time estimate > Run workflow (it first waits
for an image build still running; a newer run cancels one still going). Each study runs its full design on 4% and
6% of its datasets (at least 2 per cell; 1a and 1b in full), each on its own
4-core runner; `R/timing/estimate_run_time.R`
separates the fixed part of a study's cost from the part per dataset and reads
off the full run. The run's summary page shows the hours per study on 4, 8, 16
and 32 cores.

The same can be run anywhere, a study at a time:

```bash
VALIDATION_TIMING=0.04 Rscript R/timing/time_study.R 3a
VALIDATION_TIMING=0.06 Rscript R/timing/time_study.R 3a
Rscript R/timing/estimate_run_time.R
```

Timing runs write only to `cache/timing/`, so they never mix with a real run.

## On a cluster

The run uses one machine with many cores: the pipeline starts a worker for
every 4 cores it is given (each piece of work runs Stan's 4 chains in
parallel) and keeps every finished piece in `_targets/`. On a cluster that is
one job on one node, resubmitted until the pipeline is done. (Several nodes
would need crew's cluster controllers instead; not set up.)

**1. The project, at the commit to run.** Clone it onto storage that survives
between jobs (your home or project space, not a node's local scratch):

```bash
git clone https://github.com/dm-k1/flatcor-validation.git
cd flatcor-validation
```

**2. The image, as an Apptainer file.** Clusters seldom allow Docker, but
Apptainer (formerly Singularity) runs the same image. Pull it once; the
`.sif` file is then fixed for the whole run:

```bash
apptainer pull flatcor-validation.sif docker://ghcr.io/dm-k1/flatcor-validation:main
```

**3. A job script**, here for SLURM (`run_validation.sbatch`, next to the
project; adjust the partition, cores, memory and time to the cluster):

```bash
#!/bin/bash
#SBATCH --job-name=flatcor-validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=33
#SBATCH --mem=66G                  # 2 GB per core
#SBATCH --time=5-00:00:00
#SBATCH --output=cache/logs/slurm-%j.log

cd "$SLURM_SUBMIT_DIR"
apptainer exec --bind "$PWD":/project --pwd /project flatcor-validation.sif \
  Rscript run_all.R
```

```bash
mkdir -p cache/logs          # SLURM opens the log before the script runs
sbatch run_validation.sbatch
```

- **Cores.** The pipeline uses the job's allocation (`--cpus-per-task`), not
  all cores of the node (`parallelly::availableCores()` in `_targets.R`),
  and starts (cores - 1) / 4 workers, rounded down: one core is for the
  pipeline itself. So ask for 4 × workers + 1: 33 cores give 8 workers,
  where 32 would give 7 and leave 3 idle.
- **Memory.** Not measured yet for a full run. Bayes factor fits no longer
  keep per-observation draws, so a worker needs well under a gigabyte for
  most pieces; 2 GB per core is a safe start. After the first job,
  `sacct -j <job id> --format=MaxRSS,Elapsed` shows what it used.
- **Time.** The run-time estimate gives the core hours; divided by the cores,
  that is the wall time, longer than most clusters allow for one job. When a
  job reaches its limit, submit the same script again: the pipeline
  continues where it was. It refuses to continue with another flatcor than
  the run started with (the image fixes flatcor, so this only happens with a
  new `.sif`), and `VALIDATION_FRESH=1` (as `apptainer exec --env
  VALIDATION_FRESH=1 ...`) starts from scratch.
- **Afterwards.** The results are in `results/` (one file per study, and the
  figures). Copy them back and commit them; the notebook and the docs read
  them from there.

A SLURM array script was here until 17 November 2026 (in the git history);
one job on one node replaces it.
