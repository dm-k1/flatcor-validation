## Estimated run time of the full validation

18 November 2026: run-time estimate run
[44223161415](https://github.com/dm-k1/flatcor-validation/actions/runs/44223161415)
(commit a92e096, flatcor 07c7b55, fractions 4% and 6% of the datasets) for
1a to 3a; for 3b a calculation from timed Bayes factors (below), as its
timing runs do not fit in a GitHub runner's 6 hours. A rough estimate,
about ±30-40%: enough to choose where to run, not a number to report. The
full run's own store records the time of every step; report that one.

**Total: about 700-800 core hours** of a 2021 server core (AMD EPYC 7763,
GitHub's most common runner), of which 3b is about 600.

| Study | Core hours | How |
|---|---|---|
| 1a | 0.1 | measured in full |
| 1b | 3.7 | measured in full (EPYC 9V74) |
| 2a | 2.3 | line through 4% and 6% |
| 2b | 13 | line through 4% and 6%; likely low (the 6% run had a faster processor) |
| 3a | up to 164 | scaled up from 4% on an EPYC 7763 (98 min), which overestimates (the fixed part is scaled too) |
| 3b | about 600 | calculated, see below |
| **all** | **about 700-800** | |

### Why the range for 3a, and why "of an EPYC 7763"

GitHub's runners have several kinds of processors: this estimate ran on an
Intel Xeon 8573C and AMD EPYC 7763, 9V45 and 9V74. 3a's 4% run happened to
get the old EPYC 7763 and its 6% run a much newer EPYC 9V45, so the larger
run took less time and no line could be drawn through the two. Scaled up
from each run alone, 3a is 73 core hours of the fast processor and 164 of
the slow one. All numbers here are given for the EPYC 7763, the processor
the Mac below was compared with, so 3a counts as up to 164.

Since this run, every timing run also times a reference task
(`R/timing/reference_task.R`), which divides out the runner's speed: a
future estimate would not have this problem.

### 3b, calculated

3b at 4% is about 500 Bayes factors (each two Stan fits of 4 x 10,000
draws): 372 for the rates (3b.3-3b.6; 84 at n = 50 and 96 each at n = 200,
1000 and 5000) and 120 for growing k (3b.8-3b.9; up to 90 predictors at
n = 360 and 45 at n = 2000). One Bayes factor of each kind, timed on an
M-series Mac (one worker, 4 chains):

| Case | Seconds |
|---|---|
| n = 50, 8 vs 2 predictors | 7 |
| n = 1000, 8 vs 2 | 15 |
| n = 5000, 8 vs 2 | 58 |
| n = 5000, 3 vs 2 | 72 |
| n = 360, 90 vs 3 | 58 |
| n = 2000, 45 vs 3 | 51 |

So 3b at 4% takes about 2.5 hours for the rates and 1 hour for growing k
(about 30 s per Bayes factor on average): about 3.5 hours on the Mac. The
full 3b has 25 times the datasets: about 90 hours on the Mac. The Mac ran
study 2a at 4% in 1.0 minutes where an EPYC 7763 runner took 1.8, so it is
about 1.7 times as fast; 90 hours x 1.7 x 4 cores is about 600 core hours of
an EPYC 7763. Left out: 3b.1-3b.2 and 3b.7, which fit no models and are
cheap. Uncertain: single timings, n = 200 not timed, the growing-k average
estimated, and data where the larger model is true can cost more than null
data. The same calculation puts 3b at 4% at about 6 hours on a runner, which
is why its timing runs ran into GitHub's limit; splitting 3b's timing into
parts (rates at 2% and 4%, growing k at 4% and 6%) would measure it, but
would not change where or how the validation is run.

### What it means

| Where | Wall time | Cost |
|---|---|---|
| M-series Mac, one worker (4 of 8 cores) | about 5 days | none |
| Azure HB120rs_v3 (120 cores of the EPYC 7763's generation; 29 workers), spot | about 7 hours | about $6-8 (spot price $0.95 per hour, West Europe, 18 November 2026) |

On the cloud machine, the end of each study takes longer than the average
suggests: the last long pieces (3b at n = 5000) run while other workers
have nothing left. Splitting the run over machines, one study each
(`Rscript run_all.R 3b` on one, the rest on another), shortens that.

The numbers of datasets of 3a and 3b are still provisional ("set by cost"
in their `design.R`): chosen when the full run looked like about 1,200 core
hours. At this cost they can be set by their Monte Carlo errors instead.
