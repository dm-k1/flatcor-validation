# Running the validation elsewhere, and how long it takes

The validation is one pipeline ([`_targets.R`](../_targets.R)); `Rscript
run_all.R` runs it (only the simulation: the figures come after, from the
results files, with `Rscript make_figures.R`, on any machine). It saves every finished piece as it
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

**For your own machine.** GitHub's runners differ (several processor types),
so every timing run first times a reference task: one fixed flatcor Bayes
factor (`R/timing/reference_task.R`, about a minute). The estimate counts each
timing in multiples of its own runner's reference, so a faster or slower
runner divides out. To read the estimate off for the machine you plan to
use, run the reference task there and give its seconds to the estimate (the
artifact `run-time-estimate` holds the timing runs, in `timing/`):

```bash
Rscript R/timing/reference_task.R          # e.g. "8.9 seconds" on an M-series Mac
Rscript R/timing/estimate_run_time.R timing --reference=8.9
```

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
- **Afterwards.** The results are in `results/studies/` (one file per
  study). Copy them back, draw the figures with `Rscript make_figures.R`,
  and commit both; the notebook and the docs read them from there.

A SLURM array script was here until 17 November 2026 (in the git history);
one job on one node replaces it.

## On AWS

[`aws/user-data.sh`](aws/user-data.sh) runs the validation on one EC2 spot
instance (Amazon Linux 2023), given to it as user data with five values
filled in: the S3 bucket, a name for the run, smoke test or not, the studies,
and a safety cutoff in hours. It installs Docker, clones this repository,
runs `run_all.R` in the published image (only the simulation), copies the
results and logs to `s3://<bucket>/<run name>/`, and shuts down. The figures
are drawn afterwards from the downloaded results: `Rscript make_figures.R`. Launched with shutdown behaviour
"terminate", the instance then deletes itself.

AWS can take a spot instance back at any moment, its disk included. A real
run therefore copies the pipeline's store to `s3://<bucket>/<run name>/store`
every 15 minutes, and a new instance launched with the same run name first
copies it back and continues from there. An interruption costs at most 15
minutes of work.

[`aws/keep_running.sh`](aws/keep_running.sh), run in the background on your
own computer, launches that new instance: every 2 minutes it checks the
run's instance, and when it is gone and the run's log in S3 does not say
the run finished, it launches the fleet again (at most 12 times). When the
log says `run_all.R` stopped on an error in the code (status 1, "Execution
halted"), it stops instead and shows the end of the log, since a new
instance would stop on the same error.

The instance needs a role that may read and write the bucket. With 32 cores
the pipeline runs 7 workers.
