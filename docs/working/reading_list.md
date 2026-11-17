# Reading list

What to read before presenting or writing up this project, in order of
priority. Each entry says why it matters here and what to look for. Status:
**not read**, **partly read** (a summary or some sections), or **read**. All
citations were checked against the publisher or a library record.

What the project now claims (see [appendix B](../appendix/B_relation_to_other_priors.md),
"What is new here, and what is not"): the three priors are one family,
$R^2 \sim$ Beta($a$, 1) with a uniform direction, computed with one Stan
model and validated side by side in the same way, and each works correctly
as a Bayesian prior (the validation, studies 1a to 3b). What the results add beyond
that is **TO REFRESH** after the full run. The earlier claims about the
conditional test went with that test (removed in flatcor 2.0.0); the entries that bore
on them are kept in section 4, for the later behaviour study.

---

## 1. Read before claiming anything is new

| # | Reference | Why, and what to look for | Status |
|---|---|---|---|
| 1 | Li, Y., & Clyde, M. A. (2018). Mixtures of *g*-priors in generalized linear models. *Journal of the American Statistical Association*, *113*(524), 1828–1845. https://doi.org/10.1080/01621459.2018.1469992 | The tCCH family: whether flat R's hyperprior (a ₁F₁ in −*n*/(2*g*)) is or is not a member. The docs only say it was not found to be one. | not read |
| 2 | Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008). Mixtures of *g* priors for Bayesian variable selection. *Journal of the American Statistical Association*, *103*(481), 410–423. | The *g*-prior mixtures the three priors belong to, Zellner–Siow, and information consistency. Look for: their conditions for information consistency (compare with *n* ≥ *k* + 3 here) and model-selection consistency. | not read |
| 3 | Gabry, J., & Goodrich, B. Estimating regularized linear models with rstanarm (`stan_lm` vignette). https://mc-stan.org/rstanarm/articles/lm.html | That flat partial r is `stan_lm` with η = 1, and how they derive *R*² ~ Beta(*K*/2, η). | not read |
| 4 | Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022). Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior. *Journal of the American Statistical Association*, *117*(538), 862–874. https://doi.org/10.1080/01621459.2020.1825449 | The comparison with R2-D2 in appendix B. | not read |
| 5 | Rouder, J. N., & Morey, R. D. (2012). Default Bayes factors for model selection in regression. *Multivariate Behavioral Research*, *47*(6), 877–903. https://doi.org/10.1080/00273171.2012.734737 | The source of the JZS regression default in BayesFactor and JASP, a reference method in studies 2b and 3b. Look for: how their Bayes factor is computed (`integrate()` over *g*), and its stated precision. | not read |

## 2. Know well enough to defend the validation

| # | Reference | Why, and what to look for | Status |
|---|---|---|---|
| 6 | Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria for Bayesian model choice with application to variable selection. *Annals of Statistics*, *40*(3), 1550–1577. | The criteria checked in studies 1b and 3b, and the robust prior. Look for: the exact definitions of each criterion, so that "met" in the docs means what the paper means. | not read |
| 7 | Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788. | The calibration of the posterior in study 2a. | not read |
| 8 | Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S. (2023). Workflow techniques for the robust use of Bayes factors. *Psychological Methods*, *28*(6), 1404–1426. https://doi.org/10.1037/met0000472 | The calibration of the Bayes factor in study 3a. | not read |
| 9 | Johnson, V. E., & Rossell, D. (2010). On the use of non-local prior densities in Bayesian hypothesis tests. *Journal of the Royal Statistical Society: Series B*, *72*(2), 143–170. https://doi.org/10.1111/j.1467-9868.2009.00730.x | Why evidence for a true null grows only as a power of *n* under local priors; the frame for study 3b's rates (*a* log *n*). Look for: the rate result and its conditions. | not read |
| 10 | Mulder, J., Berger, J. O., Peña, V., & Bayarri, M. J. (2021). On the prevalence of information inconsistency in normal linear models. *TEST*, *30*(1), 103–132. | That the information-consistency statements in appendix A agree with it. | not read |
| 11 | Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. *Journal of Multivariate Analysis*, *100*(9), 1989–2001. | The LKJ distribution and its vine construction: why each partial correlation is flat under LKJ(1). | not read |
| 12 | Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology. *Journal of Mathematical Psychology*, *72*, 19–32. | Jeffreys' flat prior on a correlation, the starting point of all three priors at *k* = 1. | not read |
| 13 | Barnard, J., McCulloch, R., & Meng, X.-L. (2000). Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage. *Statistica Sinica*, *10*(4), 1281–1311. | The marginally uniform prior; the docs cite it for flat single correlations only. | not read |

## 3. Background

| # | Reference | Why | Status |
|---|---|---|---|
| 14 | Consonni, G., Fouskakis, D., Liseo, B., & Ntzoufras, I. (2018). Prior distributions for objective Bayesian analysis. *Bayesian Analysis*, *13*(2), 627–679. https://doi.org/10.1214/18-BA1103 | The survey of objective priors. Its model-comparison part does not discuss priors on *R*² or on (partial) correlations. | read (sections 3.1–3.6) |

## 4. For the later behaviour study

How the test of one coefficient behaves as predictors are added, and the
claims about the conditional test, which flatcor no longer has.

| # | Reference | Why, and what to look for | Status |
|---|---|---|---|
| 15 | Som, A., Hans, C. M., & MacEachern, S. N. (2016). A conditional Lindley paradox in Bayesian linear models. *Biometrika*, *103*(4), 993–999. | A shared scale makes the evidence for one coefficient depend on the others: the main effect the behaviour study measures. | not read |
| 16 | Pramanik, S., & Johnson, V. E. (2024). Efficient alternatives for Bayesian hypothesis tests in psychology. *Psychological Methods*, *29*(2), 243–261. https://doi.org/10.1037/met0000482 | Default tests confirm true nulls too slowly, for a psychology audience; how the number of covariates enters. | not read |
| 17 | Schönbrodt, F. D., & Wagenmakers, E.-J. (2018). Bayes factor design analysis: Planning for compelling evidence. *Psychonomic Bulletin & Review*, *25*(1), 128–142. | Operating characteristics (how often each answer occurs), the measures of the behaviour study. | not read |
| 18 | Kucharský, Š., Wagenmakers, E.-J., van den Bergh, D., & Ly, A. (2026). Analytic posterior distribution and Bayes factor for Pearson partial correlations. *TEST*, *35*(2), 535–558. | The test of one partial correlation with a flat prior (α = 1), with random predictors. Look for: how the number of control variables affects the evidence. | partly read |
| 19 | Datta, S., & Johnson, V. E. (2026). Bayes factor functions for testing partial correlation coefficients. *Journal of Mathematical Psychology*. arXiv:2503.10787. | The non-local alternative for partial correlations. | not read |
| 20 | Wetzels, R., & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test for correlations and partial correlations. *Psychonomic Bulletin & Review*, *19*(6), 1057–1064. | A JZS-based test of one partial correlation. | not read |
| 21 | Consonni, G., & Veronese, P. (2008). Compatibility of prior specifications across linear models. *Statistical Science*, *23*(3), 332–353. https://doi.org/10.1214/08-STS258 | Marginalisation, conditioning and projection for deriving priors of nested models from one full-model prior; bears on any test of one coefficient derived from the joint prior. Preprint: arXiv:1102.2981. | not read |
| 22 | Dawid, A. P., & Lauritzen, S. L. (2001). Compatible prior distributions. In E. I. George (Ed.), *Bayesian methods with applications to science, policy and official statistics* (pp. 109–118). Office for Official Publications of the European Communities. | The origin of "compatible priors". Library copy needed. | not read |
| 23 | Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010). Bayesian hypothesis testing for psychologists: A tutorial on the Savage–Dickey method. *Cognitive Psychology*, *60*(3), 158–189. | The Savage–Dickey ratio and the condition under which it equals the Bayes factor. | not read |
| 24 | Steel, M. F. J. (2020). Model averaging and its use in economics. *Journal of Economic Literature*, *58*(3), 644–719. https://doi.org/10.1257/jel.20191385 | *g*-priors and priors over model size in model averaging, if variable selection returns. | not read |

---

**Suggested order:** 6, 7 and 8 first (they define what the validation
checks), then 1, 2 and 5 (whether the family and its route are new), then the
rest of sections 1 and 2; section 4 when the behaviour study starts.
