# Appendix A. How the prior calculates

This document covers all three priors: flat partial r
(`prior = "flat_partial_r"`), flat R (`"flat_r"`) and flat R-squared
(`"flat_r_squared"`). They form one family, $R^2 \sim \mathrm{Beta}(a, 1)$ with
a uniform direction and $a = k/2$, ½ or 1, and they differ only in factor 1
below. Where a derivation is written out for one shape (flat R, $a = ½$, in
sections 7 and 8), the general case follows the same steps. Why these priors:
[`1_why_a_choice.md`](https://github.com/dm-k1/Objective-bayes/blob/main/docs/1_why_a_choice.md) and
[`2_the_three_priors.md`](../2_the_three_priors.md).

Checks are named by the study of the validation they belong to (3b.2 is study
3b, item 2; [`validation_tests.md`](../validation_tests.md) gives each study's design).
Places marked **TO REFRESH** are to be filled in from the full run.

The prior is a product of three factors:

$$
p(\beta \mid \sigma) \;=\; \underbrace{P(\text{total})}_{1}\;\times\;
\underbrace{P(\text{this combination} \mid \text{total})}_{2}\;\times\;
\underbrace{\text{unit conversion}}_{3}
$$

1. P(total) is how likely the total effect is: $R^2 \sim \mathrm{Beta}(a, 1)$,
   which is flat on the correlations ($a = k/2$), flat on $R$ ($a = ½$) or flat
   on $R^2$ ($a = 1$).
2. P(this combination | total) says that, given the total, every way of
   dividing it over the predictors is equally likely.
3. The unit conversion takes the effect from noise units to the units of the
   predictors and the outcome.

Factors 1 and 2 carry the judgements; the reasons are in
[`1_why_a_choice.md`](https://github.com/dm-k1/Objective-bayes/blob/main/docs/1_why_a_choice.md), section 3, and [`2_the_three_priors.md`](../2_the_three_priors.md), "Dividing the total". Factor 3 carries
none. It only converts units, so that the prior becomes a density per point of
$\beta$, which is what estimation software needs.

For any shape $a$, the Bayes factor of a model against the intercept-only
model is the likelihood ratio of the model's $k$ columns against none,
averaged over factor 1 (section 7, step 4). For flat R it is also the
one-predictor Bayes factor averaged over every direction the effect could take
(steps 1 to 3), and every prior can be written as a mixture of g-priors
(section 8). These derivations say what the Bayes factor is; the package does
not compute it from them. flatcor estimates each model's marginal likelihood
from its Stan model (section 6) by bridge sampling, and divides (section 7,
"How the package computes it"); the formulas serve as independent references
in flatcor's tests.

| Sections | What | Used by |
|---|---|---|
| 1–5 | the density, factor by factor | `draw_coefficients_from_prior()` |
| 6 | how estimation samples it, and the marginal likelihood | `fit_regression()`, `bayes_factor_models()` |
| 7, 9 | the Bayes factor for whole models, and how the package computes it | `bayes_factor_models()` |
| 8 | the same prior as a mixture of g-priors | the references in flatcor's tests |
| 10 | testing a term; the likelihood ratio behind the Bayes factor | `bayes_factor_models()` |
| 11–12 | the default-prior criteria; assumptions and limits | |

How this document fits into the project: [the README](../../README.md).

---

## 1. Model and notation

The model is a linear regression with $n$ observations and $k$ predictors:

$$
y = \alpha \mathbf{1} + X_c \beta + \varepsilon, \qquad \varepsilon \sim N(0, \sigma^2 I_n),
$$

where $X_c$ holds the **centred** predictor columns (each column has mean
zero), so that the intercept $\alpha$ is the mean of $y$ and separate from the
slopes.

| Symbol | Meaning |
|---|---|
| $\Sigma = X_c^\top X_c / n$ | covariance matrix of the predictors |
| $\beta^\top \Sigma \beta$ | variance explained by the predictors (this accounts for correlated predictors by itself) |
| $\sigma^2$ | unexplained variance (the noise) |
| $f^2 = \beta^\top \Sigma \beta / \sigma^2$ | explained divided by unexplained variance (**Cohen's $f^2$**) |
| $f$ | the **length of $\beta$ in noise units**: one number for the whole vector |
| $\rho^2 = f^2 / (1 + f^2)$ | the **population** share of explained variance; $R = \sqrt{\rho^2}$ |
| $R^2$ (sample) | the share of explained variance in the fitted regression |
| $\mathrm{SST} = \sum_i (y_i - \bar y)^2$ | total sum of squares |

Effects are measured in **noise units** (= with the unexplained spread
$\sigma$ as the ruler, as Cohen's $d$ divides a difference by a standard
deviation). For example, $f = 1$ means that the predictors explain as much variance as the
noise leaves unexplained, so $\rho^2 = .5$.

The prior conditions on the observed predictors, which are treated as fixed,
as in Zellner's g-prior and the JZS default. The intercept and the noise get
the usual neutral priors, as in JZS: $p(\alpha) \propto 1$ and
$p(\sigma) \propto 1/\sigma$. These priors are improper, but every model being
compared shares them, so their constants cancel in Bayes factors.

---

## 2. Factor 1: the total, R-squared ~ Beta(a, 1)

Factor 1 sets how likely each size of the total effect is. The population
share $\rho^2 = R^2$ follows $\mathrm{Beta}(a, 1)$, so $P(R^2 \le x) = x^a$,
$P(R \le u) = u^{2a}$, and the density on $R$ is $2a\,R^{2a-1}$.

| Prior | $a$ | Density on $R$ | In words |
|---|---|---|---|
| flat partial r | $k/2$ | $k R^{k-1}$ | every correlation matrix of $(y, x_1 \ldots x_k)$ equally likely |
| flat R | ½ | 1 | $R$ flat for every $k$: $P(R \le u) = u$ |
| flat R-squared | 1 | $2R$ | $R^2$ flat for every $k$ |

With one predictor, flat partial r and flat R are the same prior ($a = ½$),
Jeffreys' flat prior on the correlation.

Estimation needs this density per unit of the length $f$ instead of per unit
of $R$. With $R = f / \sqrt{1 + f^2}$:

$$
\frac{dR}{df} = (1 + f^2)^{-3/2}
\qquad\Longrightarrow\qquad
p(f) = 2a\left(\frac{f^2}{1+f^2}\right)^{a - 1/2} (1 + f^2)^{-3/2}, \qquad f > 0,
$$

and $P(|f| < f_0) = \bigl(f_0^2/(1+f_0^2)\bigr)^a$. For flat R ($a = ½$) the
density is $(1+f^2)^{-3/2}$. Equal slices of $R$ cover longer and longer
stretches of $f$, because $R$ levels off near 1, so the same probability is
spread more thinly per unit of $f$. (As a check,
$\int_0^\infty (1+f^2)^{-3/2} df = [\,f/\sqrt{1+f^2}\,]_0^\infty = 1$.) For flat
R-squared ($a = 1$) the density is $2f(1+f^2)^{-2}$, which is 0 at $f = 0$.

All three have the same tail, $p(f) \approx 2a f^{-3}$ for large $f$. This is
why they share the same boundary for information consistency (section 11).

---

## 3. Factor 2: every split of the total equally likely

Factor 2 spreads a given total evenly over all the ways it can be divided
among the predictors. Many coefficient vectors have the same length $f$: with
two predictors they form a circle, with $k$ predictors the surface of a sphere
of radius $f$ (in noise units, with the overlap between predictors taken into
account). Every combination on that sphere is equally likely, so the
probability of one combination is 1 divided by the size of the set of all
combinations:

$$
P(\text{this combination} \mid f) = \frac{1}{A_k \, f^{k-1}},
\qquad A_k = \frac{2\pi^{k/2}}{\Gamma(k/2)} .
$$

| $k$ | All combinations with length $f$ | Size |
|---|---|---|
| 1 | two points, $+f$ and $-f$ | 2 |
| 2 | a circle | $2\pi f$ |
| 3 | a sphere surface | $4\pi f^2$ |

For $k = 1$, factors 1 and 2 together give $\tfrac12 (1+f^2)^{-3/2}$ for flat
R and flat partial r. This is the flat-correlation prior on the effect-size
scale (figure 1 in `1_why_a_choice.md`, in the paper repository Objective-bayes), a $t$ distribution with 2 degrees
of freedom and scale $1/\sqrt2$ on $\beta/\sigma$. For flat R-squared they give
$|f|(1+f^2)^{-2}$.

![](../figures/4_three_factors.png)

The density has a peak at $\beta = 0$ for flat R with $k \ge 2$ and for flat
R-squared with $k \ge 3$, and none for flat partial r. For $k \ge 2$ the
size $A_k f^{k-1}$ goes to 0 as $f \to 0$, and factor 1 goes to 0 like
$f^{2a-1}$, so the density per point behaves like $f^{2a-k}$ near $\beta = 0$.
That is an infinitely high but integrable peak for flat R ($f^{1-k}$) and for
flat R-squared with $k \ge 3$ ($f^{2-k}$), and no peak for flat partial r
($f^0$). A sampler that moves through $\beta$ directly misses this peak for
larger $k$, without any warning. For that reason estimation samples $R$ and
the direction instead (section 6).

---

## 4. Factor 3: unit conversion

Factor 3 turns the density per point in noise units, which factors 1 and 2
give, into a density per point in the units of the data, which is where the
coefficients $\beta$ live. Changing the axes stretches or squeezes space, and
the density changes by the inverse of that stretch:

$$
\text{unit conversion} = \frac{\sqrt{\det \Sigma}}{\sigma^k} .
$$

The conversion is the same for every $\beta$. In the units of the data, the
circles of factor 2 become tilted ellipses, because predictors differ in
spread and are correlated.

---

## 5. The whole density

Putting the three factors together, on the log scale:

$$
\log p(\beta \mid \sigma)
= \underbrace{\log(2a) + \left(a - \tfrac12\right)\log\frac{f^2}{1+f^2} - \tfrac32 \log(1 + f^2)}_{1}
\;\underbrace{-\,\log A_k - (k-1)\log f}_{2}
\;+\underbrace{\tfrac12 \log\det\Sigma - k \log \sigma}_{3},
\qquad f^2 = \frac{\beta^\top\Sigma\beta}{\sigma^2} .
$$

Drawing from the prior follows the same three factors in order.
`draw_coefficients_from_prior()` draws $R = U^{1/(2a)}$ with $U$ uniform
(factor 1), a direction with no preference in explained variance (factor 2),
and scales it to length $f$ in the units of the data (factor 3). Counting those
draws in small squares reproduces the density of the formula (figure below;
the predictors are correlated, so factor 3 is tested too).

![](../figures/5_formula_vs_draws.png)

---

## 6. Estimation: sampling in the prior's own coordinates

The density in section 5 describes the prior, but sampling it directly in
$\beta$ works poorly. Near $\beta = 0$ the density has the peak described in
section 3, and a sampler that moves through $\beta$ steps over it without
warning; with 8 predictors it draws far too few small effects.

[`inst/stan/flatcor.stan`](https://github.com/dm-k1/flatcor/blob/main/inst/stan/flatcor.stan) in flatcor therefore samples the parameters from which the prior
is built:

| Parameter | Prior | Factor |
|---|---|---|
| `R_magnitude` $= \lvert R\rvert$ | density $2a\,R^{2a-1}$ on (0, 1) | 1 |
| `direction` $= z/\lVert z\rVert$ | uniform on the sphere | 2 |
| `sigma`, `alpha` | $1/\sigma$, flat | |

With several predictors the sign of each coefficient lives in the direction.
With one predictor the sign of $R$ is + or − with probability ½: the model
averages the likelihood over both signs and draws the sign afterwards, from
the probability the likelihood gives it (flatcor 1.0.1). Sampling a signed $R$
on (−1, 1) instead, as flatcor 1.0.0 did, puts a zero of the density at
$R = 0$ whenever $a > ½$: under flat R-squared the density is $\lvert R\rvert$,
and the sampler has to cross that zero at every change of sign, which gives
divergent transitions.

It then computes $\beta = \sigma f\, L^{-\top} d$, with $f = R/\sqrt{1-R^2}$,
$d$ the direction and $L$ a square root of $\Sigma$ (here the Cholesky factor).
This gives $\beta^\top \Sigma \beta = \sigma^2 f^2$, so the length is correct
(factor 3). A matrix root is needed at this point, because it turns a uniform
direction in explained variance into predictor units; which root is used does
not matter.

The direction is $z/\lVert z\rVert$ for a free vector $z$. Any density for $z$
that depends only on its length makes the direction uniform, and the model
uses one that keeps $\lVert z\rVert$ near 1. Stan's built-in `unit_vector`
does the same with a standard normal $z$, which lets $\lVert z\rVert$ come
close to 0; there the direction is undefined and the sampler can struggle.

The same Stan model serves all three priors, because the prior enters only
through factor 1. The data value `r_squared_shape` is $a$, and one line adds
`target += log(2 * r_squared_shape) + (2 * r_squared_shape - 1) * log(R_magnitude)`.
Factors 2 and 3 are the same for all three. Every density is written with all
its normalising constants (including that of $z$'s density, which depends on
$k$), so that the model's marginal likelihood can be estimated from it
(section 7). With no predictors ($K = 0$) it is the intercept-only model.

In these coordinates the prior has no peak. Study 1a (`studies/1a_prior_sampling/`)
checks that the prior is reproduced for $k$ = 1, 3 and 8 for each prior, and
study 2a that the posterior is calibrated (simulation-based calibration, all
three priors).

---

## 7. The Bayes factor for whole models

This section derives the Bayes factor of the model with $k$ predictors against
the intercept-only model,
$\mathrm{BF}_k = p(y \mid k \text{ predictors}) / p(y \mid \text{intercept only})$.
Steps 1 to 3 derive it for flat R ($a = ½$) through g-priors; section 8
extends that derivation to every $a \le k/2$. Step 4 gives the same Bayes
factor for every $a$ as one integral over the effect size. The package
estimates it by bridge sampling instead ("How the package computes it",
below); these formulas are its independent references.

The idea of steps 1 to 3 is to fix a direction first. Along one direction the model is a
regression on a single constructed predictor, whose Bayes factor is known.
Averaging over directions then gives the Bayes factor for $k$ predictors.

### Step 1: one predictor along a fixed direction

This step rewrites the model so that a fixed direction becomes a single
predictor. Work with the **fitted values** $X_c\beta$ instead of the
coefficients. Take any orthonormal basis $Q$ ($n \times k$) of the columns of
$X_c$, for example `qr.Q(qr(X_c))`, the QR decomposition that `lm()` uses.
Every vector of fitted values is $X_c\beta = Q w$ for some $w$ with $k$
entries, and

$$
\lVert w \rVert^2 = \lVert X_c \beta \rVert^2 = n\, \beta^\top \Sigma \beta .
$$

In $w$, the explained variance is therefore an ordinary squared length, and
the flat R prior reads: length $\lVert w \rVert = \sigma\sqrt n\, f$, and
direction $w / \lVert w \rVert$ uniform on the sphere. The choice of
orthonormal basis does not matter, because a uniform direction looks the same
in every orientation.

For a fixed direction $d$ (a unit vector), the signal is
$\sigma \sqrt n\, f\,(Qd)$, that is, one predictor $v = Qd$ with
$\lVert v \rVert = 1$. If $f$ carries a sign ($d$ and $-d$ are equally likely),
the coefficient in noise units has density $\tfrac12(1 + f^2)^{-3/2}$, the
one-predictor flat-correlation prior.

The share of explained variance in this one-predictor regression is

$$
R_d^2 = \frac{(v^\top y_c)^2}{\mathrm{SST}} = R^2\, c^2 ,
$$

where $c$ is the **alignment** (= the cosine of the angle between $d$ and the
direction of the data, the least-squares fit $Q^\top y_c$). Given
$d$, the data enter only through $R^2$ and $c$.

### Step 2: the one-predictor Bayes factor, through g

The one-predictor Bayes factor needs an integral over $\sigma$ (and $\alpha$),
and that integral has a known closed form when the prior on the coefficient
is normal. This is Zellner's g-prior, in which $g$ is the prior variance of
the coefficient relative to the noise. The flat-correlation prior is therefore
written as a mixture of normal priors:

$$
f \mid g \sim N(0, g/n), \qquad g \sim \text{Inverse-Gamma}(1, \; n/2),
$$

since $\int (g/n)^{-1/2} e^{-nf^2/(2g)}\, g^{-2} e^{-n/(2g)}\, dg \propto (1+f^2)^{-3/2}$.
For each $g$ the g-prior Bayes factor is known (Liang et al., 2008), with
$R_d^2$ in place of $R^2$:

$$
\mathrm{BF}(g, d) = (1+g)^{(n-2)/2}\,\bigl(1 + g(1 - R_d^2)\bigr)^{-(n-1)/2}
= (1+g)^{-1/2} \left(1 - R^2 c^2 \frac{g}{1+g}\right)^{-(n-1)/2} .
$$

The Zellner–Siow (JZS) default uses Inverse-Gamma(½, $n s^2/2$) instead,
which is a Cauchy prior on the effect with scale $s$ (0.354 in BayesFactor and
JASP, 1 for Zellner–Siow). The shape of this hyperprior, 1 instead
of ½, is what separates "flat on the correlation" from "Cauchy on the effect
size".

### Step 3: average over g and over directions

Averaging the Bayes factor of step 2 over $g$ and over directions gives the
Bayes factor for $k$ predictors. For a uniform direction in $k$ dimensions,
the alignment $c$ with any fixed direction has density

$$
p_k(c) \propto (1 - c^2)^{(k-3)/2}, \qquad -1 < c < 1 ;
$$

its square $c^2$ is the **share** of [`2_the_three_priors.md`](../2_the_three_priors.md),
"The share per predictor", which follows Beta(½, (k−1)/2). Averaging gives

$$
\boxed{\;
\mathrm{BF}_k(R^2, n) = \int_0^\infty \int_{-1}^{1}
\mathrm{IG}(g; 1, n/2)\,(1+g)^{-1/2}
\left(1 - R^2 c^2 \frac{g}{1+g}\right)^{-(n-1)/2} p_k(c)\, dc\, dg
\;}
$$

The data enter only through $R^2$, $n$ and $k$. For $k = 1$ there are only two
directions, both with $c^2 = 1$, and the formula reduces to the one-predictor
flat-correlation Bayes factor.

### Step 4: the same Bayes factor for every shape a

The whole-model Bayes factor is also the likelihood ratio of the model's $k$
columns against none, for a length $f$ of their effect in noise units with a
uniform direction (section 10, with $q = k$, so $m = n - 2$),

$$
\frac{p(y \mid f)}{p(y \mid 0)} = e^{-nf^2/2}\, {}_1F_1\!\left(\tfrac{n-1}{2}; \tfrac k2; \tfrac{n R^2 f^2}{2}\right),
$$

averaged over factor 1. Since $R^2 = f^2/(1+f^2) \sim$ Beta($a$, 1), the
length has density $p_a(f) = 2a\, f^{2a-1}/(1+f^2)^{a+1}$, and

$$
\boxed{\;\mathrm{BF}_k = \int_0^\infty e^{-nf^2/2}\, {}_1F_1\!\left(\tfrac{n-1}{2}; \tfrac k2; \tfrac{n R^2 f^2}{2}\right) p_a(f)\, df\;}
$$

Only $p_a$ differs between the priors. At $n = k + 1$ every model fits the
data exactly ($R^2 = 1$), and the data say nothing about the effect; the
Bayes factor is then 1 for any proper prior on it (predictive matching,
Bayarri et al., 2012).

### How the package computes it

The package does not evaluate these formulas. `bayes_factor_models()` fits
each of the two models with the Stan model of section 6 (the chosen prior,
the same improper priors on the intercept and the noise in both), estimates
each model's marginal likelihood from its posterior draws by bridge sampling
(bridgesampling; Gronau, Singmann & Wagenmakers, 2020), and divides. A model
with no predictors is the intercept-only model, so the Bayes factor against it
is the same call. The result is an estimate with a Monte Carlo error, which the
package reports (bridge sampling's approximate error of the log Bayes factor).
Following the bridgesampling authors, each model gets ten times the draws used
for estimation. The formulas of this section and section 8 are independent
references for it in flatcor's tests; validation study 3a checks it end to end
by calibration.

---

## 8. The same prior as a mixture of g-priors

The family $R^2 \sim$ Beta($a$, 1) with a uniform direction is a family of
mixtures of Zellner g-priors (for $a \le k/2$). Within it, flat R is the only
member whose implied $R^2$ is the same for every $k$ *and* equals Jeffreys'
flat correlation at $k = 1$.

Zellner's g-prior is $\beta \mid g, \sigma \sim N\bigl(0,\ g\,\sigma^2 (X_c^\top X_c)^{-1}\bigr)$:
a uniform direction (in explained variance) with $f^2 = (g/n)\,\chi^2_k$. A
prior on $g$ gives a **mixture of g-priors** (Liang et al., 2008); the
Zellner–Siow (JZS), hyper-g, hyper-g/n, beta-prime, robust and intrinsic priors
are all members. Flat R is the member with

$$
\frac{g}{n} = B \cdot W, \qquad B \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
\qquad W \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right).
$$

- $W$ is the one-predictor part: for $k = 1$, $B = 1$, and this is the
  Inverse-Gamma(1, $n/2$) of section 7, step 2.
- $B$ is again the share: the squared alignment $c^2$ of section 7, and the
  share of [`2_the_three_priors.md`](../2_the_three_priors.md), "The share per predictor". Averaging over
  directions in section 7 is the same as averaging over $B$ here.
- The construction works because $f^2 = B W \chi^2_k$, and $B\,\chi^2_k$ has the
  distribution of $\chi^2_1$ (a Beta(½, (k−1)/2) share of a $\chi^2_k$ is a
  $\chi^2_1$). So $f^2 = W\chi^2_1$ for every $k$, which is exactly the
  one-predictor prior on the total.
- Being constant in $k$ does not single this prior out by itself. Any $W$
  gives an implied $R^2$ that does not depend on $k$: with
  $W \sim$ Inverse-Gamma($b$, ½) it is $R^2 \sim$ Beta(½, $b$) for every $k$
  ($W\chi^2_1$ is then a ratio of a Gamma(½) and a Gamma($b$) variable). The second condition fixes $b = 1$: at $k = 1$ the
  prior must be flat on $r$, that is, $r^2 \sim$ Beta(½, 1). More generally,
  for a g-mixture $f^2 = (g/n)\,\chi^2_k$, so the distribution of $f^2$ at each
  $k$ determines that of $g/n$ (Mellin transforms factor). Fixing the $k = 1$
  distribution to Jeffreys' and requiring it to stay the same for every $k$
  leaves exactly one hyperprior per $k$, this one.

The hyperprior density is $p(g/n = t) = \dfrac{1}{2kt^2}\, {}_1F_1\!\left(\tfrac32;
\tfrac k2 + 1; -\tfrac{1}{2t}\right)$, where ${}_1F_1$ is Kummer's confluent
hypergeometric function. That factor is why the prior is not a special case of
the tCCH family of Li and Clyde (2018), which covers hyper-g, beta-prime,
robust and intrinsic priors: their densities are products of powers and an
exponential.

The Bayes factor is then the standard g-prior Bayes factor averaged over this
hyperprior, a single integral over $g$:

$$
\mathrm{BF}_k = \int_0^\infty (1+g)^{(n-1-k)/2}\,\bigl(1+g(1-R^2)\bigr)^{-(n-1)/2}\, p(g)\, dg ,
$$

flatcor's tests evaluate it with `integrate()` as an independent reference for
the package's Bayes factor (section 7, "How the package computes it").

The table places the three priors among the other mixtures by their median
implied $R^2$ before any data are seen ($n = 100$):

| Prior | $k = 1$ | $k = 3$ | $k = 10$ | $k = 30$ |
|---|---|---|---|---|
| Zellner–Siow (scale 1) | .50 | .84 | .95 | .98 |
| hyper-g/n ($a = 3$) | .54 | .87 | .96 | .99 |
| robust (Bayarri et al., 2012) | .52 | .72 | .77 | .78 |
| intrinsic | .39 | .59 | .65 | .66 |
| **flat R** | .25 | .25 | .25 | .25 |
| **flat R-squared** ($a = 1$) | .50 | .50 | .50 | .50 |
| **flat partial r** ($a = k/2$) | .25 | .63 | .87 | .95 |

The robust and intrinsic priors also scale $g$ with $k$, and drift less than
Zellner–Siow, but their implied $R^2$ still changes with $k$.

The whole family follows from the same construction. With
$R^2 \sim$ Beta($a$, 1) and a uniform direction, the prior is the g-mixture
with $B \sim$ Beta($a$, $k/2 - a$) and the same $W$, for $a \le k/2$ (the
reference in flatcor's tests):

| Prior | $B$ | A g-mixture? |
|---|---|---|
| flat R ($a = ½$) | Beta(½, (k−1)/2) | for every $k$ |
| flat R-squared ($a = 1$) | Beta(1, $k/2 - 1$) | for $k \ge 3$; at $k = 2$, $B \equiv 1$ and it is flat partial r; at $k = 1$ no ($a > k/2$: it would need $B > 1$); there the integral over the effect itself (section 7, step 4) is the reference |
| flat partial r ($a = k/2$) | $B \equiv 1$: $g/n \sim$ Inverse-Gamma(1, ½) for every $k$ | for every $k$ |

At $a = k/2$ the prior is the one obtained when every correlation matrix of
$(y, x_1 \ldots x_k)$ is equally likely, and it is rstanarm `stan_lm`'s $R^2$
prior with η = 1. Zellner–Siow (scale 1) is exactly $R^2 \sim$ Beta($k/2$, ½)
(its row in the table above). How the three compare:
[appendix B](https://github.com/dm-k1/Objective-bayes/blob/main/docs/appendix/B_relation_to_other_priors.md).

![](../figures/6_g_prior_family.png)

---

## 9. Comparing any two models

Any two models can be compared through their Bayes factors against the same
intercept-only model, which uses the same priors on $\alpha$ and $\sigma$. For
two models A and B fitted to the same outcome and rows, nested or not:

$$
\mathrm{BF}_{AB} = \frac{\mathrm{BF}_{k_A}(R^2_A, n)}{\mathrm{BF}_{k_B}(R^2_B, n)} .
$$

The evidence for one coefficient, or for a block of columns such as a factor,
is this ratio for the models with and without it. Each model uses the prior
with the covariance of its own predictors and its own $k$.

In the package, `bayes_factor_models()` estimates each model's marginal
likelihood by bridge sampling and divides (section 7, "How the package
computes it"); the intercept-only model's marginal likelihood cancels, so it
is never needed for comparing two models with predictors.

---

## 10. Testing a term, and the likelihood ratio behind the Bayes factor

A term (one coefficient, or a factor's or an interaction's columns together)
is tested by comparing the models with and without it (section 9), each with
its own prior. This test judges the term by how it changes the total, so it
depends on how much the other predictors explain. For example, take two
datasets in which *x*₁ has the same partial correlation, .2 ($n = 100$,
$k = 3$). When the other predictors explain 5%, the flat R Bayes factor for
*x*₁ is 2.25; when they explain 60%, it is 0.70, with the same evidence about
*x*₁. This is the **conditional Lindley paradox** of Som, Hans and
MacEachern (2016). Every prior with one scale shared by all coefficients has
it, including the JZS default. Under flat R and flat R-squared the test also
depends on how many other predictors there are, because the prior for one
coefficient narrows as $k$ grows (the motivating simulation,
[`1_why_a_choice.md`](https://github.com/dm-k1/Objective-bayes/blob/main/docs/1_why_a_choice.md), section 3): these priors expect
each coefficient to take a small share of the total, on average $1/(3k)$ and
$1/(2k)$ of the variance. How much both matter in practice is for the later
study of how the priors behave in use.

### The likelihood ratio of an effect's size

Every Bayes factor here is one likelihood ratio averaged over the prior,
which is why it depends on the data only through $R^2$, $n$ and $k$. Take $q$ columns of a model with $k$ predictor columns ($q = k$: the
whole model; $q = 1$: one coefficient). Let $f$ be the length of their effect
in noise units (its $f^2$ is the variance they add beyond the other
predictors, divided by $\sigma^2$), with its direction within the $q$ columns
uniform. The other coefficients and the intercept get flat priors and the
noise gets $1/\sigma$, the same under both hypotheses. Integrating them out
gives

$$
\frac{p(y \mid f)}{p(y \mid 0)} = e^{-nf^2/2}\,\frac{\mathrm E_c\, I(\sqrt n\, r f c)}{I(0)}
= e^{-nf^2/2}\; {}_1F_1\!\left(\tfrac{m+1}{2}; \tfrac q2; \tfrac{n r^2 f^2}{2}\right),
$$

$$
I(b) = \int_0^\infty v^{m} e^{-v^2/2 + bv}\,dv, \qquad m = n - k + q - 2,
$$

with $r^2$ the columns' sample partial $R^2$ (from the models with and
without them) and $c$ the cosine between the direction of the effect and
that of the data. The ratio depends only on $r^2$, $n$, $k$ and $q$. The
average over directions is in closed form,

$$
\mathrm E_c\, e^{tc} = {}_0F_1\!\left(; \tfrac q2; \tfrac{t^2}{4}\right)
= \Gamma\!\left(\tfrac q2\right) \left(\tfrac t2\right)^{1 - q/2} I_{q/2-1}(t),
$$

a Bessel function ($\cosh t$ for $q = 1$, the two signs).

With $q = k$, averaged over the prior on $f$, this is the whole-model Bayes
factor (section 7, step 4).

Study 3a checks that each prior's Bayes factor for a term is calibrated on
data drawn from that prior: among datasets with a given BF₁₀, the share
without the term is 1/(1 + BF₁₀).

---

## 11. The criteria for a default prior (Bayarri et al., 2012)

Bayarri, Berger, Forte and García-Donato (2012) list seven criteria that a
default prior for model choice should meet. The table gives the status for
each of the three priors. "Checked" means a numerical check in study 3b
(`studies/3b_default_criteria/`; invariance in study 1b), with flatcor's
Bayes factor by bridge sampling, in the settings it covers (growing $k$ up to
90); none of the criteria is proven analytically here. **TO REFRESH:** the
results in the table, from the full run.

| # | Criterion | Meaning | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|---|
| 1 | Basic | the prior on the tested coefficients is proper, and not arbitrarily vague | **met**, by construction | **met** | **met** |
| 2 | Model selection consistency | the true model's posterior probability goes to 1 as $n$ grows | checked (3b.3–3b.6, 3b.8–3b.9): expected null rate $\tfrac q2 \log n$ in every case. **TO REFRESH:** slopes | checked (same): expected $\tfrac12 \log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise. **TO REFRESH:** slopes | checked (same): expected $\log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise. **TO REFRESH:** slopes |
| 3 | Information consistency | BF → ∞ as $R^2 \to 1$ | **met for $n \ge k + 3$** (derived below); checked (3b.2) | **met for $n \ge k + 3$**; checked (3b.2) | **met for $n \ge k + 3$**; checked (3b.2) |
| 4 | Intrinsic consistency | the prior converges to a fixed proper prior as $n$ grows | checked (3b.7) | checked (3b.7) | checked (3b.7) |
| 5 | Predictive matching | on minimal samples ($n = k + 1$) BF = 1 with the null and between models of equal size | **met** (section 7); checked (3b.1) | **met**; checked (3b.1) | **met**; checked (3b.1) |
| 6 | Measurement invariance | units of the data do not change the answer | **met**, by construction; checked (1b) | **met**; checked (1b) | **met**; checked (1b) |
| 7 | Group invariance | the common parameters get the right-Haar prior, and the prior respects the invariance | **met** | **met** | **met** |

Dimensional matching follows from predictive matching for all three priors:
the Bayes factor depends only on $R^2$, $n$ and $k$, and $R^2 = 1$ for every
model at $n = k + 1$. Group invariance holds because
$p(\alpha, \sigma) \propto 1/\sigma$ and $\beta$ is scaled by $\sigma$ in all
three.

The boundary $n \ge k + 3$ for information consistency follows from the tail
of the hyperprior. As $R^2 \to 1$, the g-prior Bayes factor of section 8 tends
to $\int (1+g)^{(n-1-k)/2}\,p(g)\,dg$, so the Bayes factor grows without limit
exactly when that integral diverges (Liang et al., 2008). For large $g$, every
member of the family has $p(g/n = t) \asymp t^{-2}$ (for flat R the confluent
hypergeometric factor tends to 1, so $p(t) \approx 1/(2kt^2)$; for
$B \sim$ Beta($a$, $k/2-a$) the bounds $0 < B \le 1$ give the same order). The
integrand therefore behaves like $g^{(n-1-k)/2 - 2}$, which diverges if and
only if $(n-1-k)/2 - 2 \ge -1$, that is, $n \ge k + 3$. On the effect scale
this is the same tail, $p(f) \approx 2a f^{-3}$ for every $a$ (section 2).
The boundary itself is not checked numerically: at $n = k + 2$ the evidence
levels off only at $R^2$ within $10^{-8}$ of 1, beyond what a sampler can
resolve.

Criterion 2 is checked for nested models only; non-nested comparisons are
listed as open in the [README](../../README.md#open-questions).

---

## 12. Assumptions and limits

- Normal errors, fixed effects and an intercept.
- The prior uses the **sample** covariance of the predictors (it conditions
  on $X$), like the g-prior and JZS.
- The improper priors on $\alpha$ and $\sigma$ are valid only for comparing
  models that share them, which all models here do.
- Not covered: group-level effects, other likelihoods and distributional
  models.
