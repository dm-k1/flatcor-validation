# How these results were made

The record of the real runs behind the files in this folder: which run made
each file, with what software, on which machines, and what happened on the
way. What the studies test and why is in
[`docs/validation_tests.md`](../docs/validation_tests.md).

## The files

| File | Made by | Commit | flatcor | Smoke test |
|---|---|---|---|---|
| `studies/study_1a.rds` | run 1a–2b | 5e91c5e | 07c7b55 | no |
| `studies/study_2a.rds` | run 1a–2b | 5e91c5e | 07c7b55 | no |
| `studies/study_2b.rds` | run 1a–2b | 5e91c5e | 07c7b55 | no |
| `studies/study_1b.rds` | run 1b, 3a, 3b | 467b18d | 07c7b55 | no |
| `studies/study_3a.rds` | run 1b, 3a, 3b | 467b18d | 07c7b55 | no |
| `studies/study_3b.rds` | run 1b, 3a, 3b | 467b18d | 07c7b55 | no |
| `figures/2a_sbc_ecdf_difference.png`, `figures/2b_bias.png`, `figures/2b_coverage.png`, `figures/2b_prediction.png` | run 1a–2b | 5e91c5e | | |
| `figures/3a_bayes_factor_calibration.png`, `figures/3b_growing_k.png` | run 1b, 3a, 3b | 467b18d | | |

Each results file records its own commit, flatcor commit, R and CmdStan
versions and finishing time (`readRDS(file)$setup`). The code of 1a, 2a and 2b
is the same at both commits: 467b18d changed only 1b, 3a, 3b and the rules
for how many datasets (below).

## Software

Every machine ran the same container image,
`ghcr.io/dm-k1/flatcor-validation:main`, built by `.github/workflows/image.yml`
from 13c6e41 (18 November 2026). Its environment, `renv.lock` and
`cluster/Dockerfile`, is unchanged from 1549f1a on, so the commits above
change the validation's code only, never its software.

| | Version |
|---|---|
| R | 4.6.1 |
| CmdStan | 2.40.0 (cmdstanr 0.9.0) |
| flatcor | 2.0.0.9000, commit 07c7b55 |
| bridgesampling | 1.2-1 |
| BayesFactor | 0.9.12-4.8 |
| SBC | 0.5.0.0 |
| targets / crew / mirai | 1.12.0 / 1.3.3 / 2.7.4 |
| posterior | 1.7.0 |
| dplyr | 1.2.1 |
| everything else | `renv.lock` |

## The runs

Both ran on AWS EC2 (Amazon Linux 2023), one machine at a time, with
`Rscript run_all.R <studies>` inside the image (`cluster/aws/user-data.sh`).
Every piece of work runs Stan's 4 chains; a machine runs (cores − 1) / 4
pieces at once. The store of finished pieces was copied to S3 every 15
minutes, and a machine taken back by AWS (spot) was replaced by one that
continued from that copy; only the pieces in progress were computed again.
Each piece draws its data and its Stan seed from its own ID, so a piece
computed again gives the same result.

**Run 1a–2b** (studies 1a, 1b, 2a, 2b; commit 5e91c5e), 19 November 2026,
00:45–04:46 CET. On an on-demand c7a.4xlarge in Stockholm (3 workers), moved
halfway to a spot c7a.8xlarge in N. Virginia (7 workers). 56 core hours. Its
1b was later replaced by the 1b of the next run.

**Run 1b, 3a, 3b** (commit 5e91c5e for round 1 of 3a and 3b, then 467b18d),
18 November evening to 19 November 17:42 CET. Spot machines in Stockholm and
N. Virginia: AMD EPYC (Zen 4: c7a and m7a, 8xlarge and 16xlarge, 7 or 15
workers) and, for part of round 1, Intel Xeon (c6i.8xlarge). AWS took a
machine back eight times; each time the run continued from the copy in S3.

- Round 1 (3a and 3b) finished at 05:25 CET. The precision rules of the time
  then asked for about 1.4 million more datasets, mostly to estimate slopes
  whose sign was already certain, and the budget was only checked between
  rounds. The run was stopped at 07:25 CET.
- The rules were revised in 467b18d (`docs/validation_tests.md`, "How many
  datasets"): Monte Carlo errors measured across datasets, a relative rule
  for direction checks, 3a decided by the data-averaged posterior, and a run
  that stops before a round it cannot afford. Every finished piece stayed
  valid; the round-2 plans of the old rules were removed.
- The run continued at 467b18d with 1b, 3a and 3b from 09:28 CET: 1b from
  scratch (so that it gets the measured errors), 3a and 3b from their
  finished pieces. After round 2 every cell had reached its targets.
- 689 core hours in this run's store, within the budget of 2,500. The count
  is a piece's duration × 4 cores, and on overloaded machines it runs high.

## How much

| Study | Datasets | Stan fits (flatcor) | Other Bayesian computations |
|---|---|---|---|
| 1a | — (the prior only) | 9 | — |
| 1b | 150 | 1,500 | — |
| 2a | 786 | 786 | — |
| 2b | 13,700 | 41,100 | 27,400 (JZS, Zellner–Siow) |
| 3a part A | 1,250 | 2,500 | — |
| 3a part B | 560 (5 runs each) | 5,600 | — |
| 3b.1, 3b.2 | 21 | 42 | — |
| 3b rates | 5,170 | 10,340 | — |
| 3b growing *k* | 11,040 | 18,080 | about 4,000 (JZS, Zellner–Siow) |
| **total** | **about 32,700** | **about 80,000** | **about 31,000** |

A Bayes factor is two Stan fits and bridge sampling: about 21,000 in all.

## Reproducing

The same commit on the same family of processor gives identical tables (in a
test of a multi-machine setup, a smoke run on one machine and on several gave
the same 32 tables). Across processor families, Stan's floating-point
arithmetic can differ in the last digits, and a sampler carries that into
different draws. An earlier run of 1a–2b in a GitHub codespace (AMD Zen 3,
commit b75a082, not part of these results) gave an identical 1a, and for 1b,
2a and 2b the same conclusions, with numbers that differ within their Monte
Carlo error.
