# Per-dataset results of the final validation

One row per dataset (and prior or method) of every study: the numbers the
summaries in `../studies/` were made from. Exported once on 20 Nov 2026 by
`export_v1.R` from the two final targets stores, without running any step
and without computing any number. The checks at the end of that script
confirm that the stored rows are the ones the results files summarise
(2b and 3a row for row; the medians of 3b's growing k equal the reported
ones). Which commit and flatcor made each study: `../RUNS.md`.

The design columns name a dataset: with the replication they are its name,
from which its data and seeds follow (`docs/validation_tests.md`,
"Randomness"). `divergent` is the number of divergent transitions;
`error` is bridge sampling's error of the log Bayes factor.

| File | Rows | One row is | Columns beyond the design |
|---|---|---|---|
| `study_1a_prior_draws.csv.gz` | 180,000 | one draw of the prior: `prior`, `design`, `quantity`, `chain`, `iteration` | `value` |
| `study_1b_changes.csv` | 450 | one dataset's change in a log Bayes factor when the data are rescaled or recoded | `expected`, `estimate`, `mcse`, `distance_in_mcse`, `bridge_error` |
| `study_2a_sbc_fits.csv` | 2,358 | one parameter of one simulation-based calibration fit (`sim_id`) | `simulated_value`, `rank`, posterior summaries, `rhat`, `ess_bulk`, `ess_tail` |
| `study_2b_fits.csv` | 82,200 | one method's fit of one dataset | `estimate`, `lower`, `upper` (interval), `prediction_error`, `divergent` |
| `study_3a_calibration.csv` | 1,275 | one dataset's log Bayes factor, data from either model (`from_model_without_x1`) | `log_BF10`, `error`, `divergent` |
| `study_3a_repeated_runs.csv` | 2,800 | one of several runs on the same dataset (`run`) | `log_BF10`, `error`, `divergent` |
| `study_3b_formula.csv` | 21 | one prior at one small *n* (items 3b.1, 3b.2) | `log_bf`, `error`, `divergent` |
| `study_3b_rates.csv` | 5,170 | one dataset's log Bayes factor (items 3b.3–3b.6) | `log_bf`, `error`, `divergent` |
| `study_3b_growing_k.csv` | 11,100 | one method on one dataset as *k* grows with *n* (items 3b.8, 3b.9) | `log_bf`, `error`, `divergent` |

Not here, as the runs did not keep them: the draws behind 3b.7 (intrinsic
consistency; only its test is in `../studies/study_3b.rds`), and the two log
Bayes factors per dataset behind 1b's changes (only the changes).
