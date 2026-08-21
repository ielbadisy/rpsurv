# rpsurv

<!-- badges: start -->
[![R-CMD-check](https://github.com/ielbadisy/rpsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ielbadisy/rpsurv/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`rpsurv` fits the Royston and Parmar (2002) family of flexible parametric
survival models: the baseline distribution is a restricted cubic spline in
log time on a chosen transformed scale (proportional hazards, proportional
odds, or probit), which gives a smooth, fully parametric, extrapolable fit
without committing to a Weibull, log-normal, or other single distribution.

The log-likelihood and its analytic gradient are evaluated in C++ via
`RcppParallel`, so fitting stays fast at large sample sizes (see the
Benchmark section below). Output mirrors `coxph()`: `print()`/`summary()`
separate interpretable covariate effects from the baseline spline nuisance
parameters, and there are `predict()`/`plot()`/residual-diagnostic methods.

All results below are real output from `rstpm2::brcancer` (686 rows, 299
events).

## Installation

```r
# install.packages("rpsurv")  # once on CRAN
# development version:
# remotes::install_github("ielbadisy/rpsurv")
```

## Fitting a model

```r
library(rpsurv)
library(survival)

data(brcancer, package = "rstpm2")
brcancer$hormon <- as.numeric(brcancer$hormon)

fit <- rpsurv(Surv(rectime, censrec) ~ hormon, data = brcancer, df = 4, scale = "hazard")
fit
#> Royston-Parmar flexible parametric survival model
#> Call:
#> rpsurv(formula = Surv(rectime, censrec) ~ hormon, data = brcancer,
#>     df = 4, scale = "hazard")
#>
#> Scale: hazard   Baseline df: 4
#> n = 686 , number of events = 299
#> Log-likelihood = -2606   AIC = 5225

summary(fit)
#> Royston-Parmar flexible parametric survival model
#> Call:
#> rpsurv(formula = Surv(rectime, censrec) ~ hormon, data = brcancer,
#>     df = 4, scale = "hazard")
#>
#> Covariate effects (hazard scale):
#>           coef exp(coef) [HR] se(coef)      z Pr(>|z|)
#> hormon -0.3641         0.6948   0.1249 -2.914  0.00356 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline spline terms:
#>                   coef  exp(coef)   se(coef)      z Pr(>|z|)
#> (Intercept) -1.960e+01  3.086e-09  3.257e+00 -6.016 1.79e-09
#> s(logt)1     2.993e+00  1.994e+01  5.966e-01  5.016 5.26e-07
#> s(logt)2    -4.429e-01  6.422e-01  5.711e-01 -0.776   0.4380
#> s(logt)3     1.357e+00  3.884e+00  8.188e-01  1.657   0.0974
#> s(logt)4    -8.841e-01  4.131e-01  4.176e-01 -2.117   0.0342
#>
#> n = 686 , number of events = 299 , parameters = 6
#> Log-likelihood = -2606   AIC = 5225   BIC = 5247
```

The summary keeps the interpretable covariate hazard ratio (`hormon`)
separate from the baseline spline coefficients, which describe the shape
of the baseline hazard and are not individually interpretable, matching
`rstpm2`/`flexsurv` convention.

## Prediction and plotting

```r
predict(fit, newdata = data.frame(hormon = c(0, 1)), times = c(500, 1000, 2000), type = "survival")
#>   id time   est
#> 1  1  500 0.829
#> 2  1 1000 0.628
#> 3  1 2000 0.408
#> 4  2  500 0.878
#> 5  2 1000 0.724
#> 6  2 2000 0.536

plot(fit, newdata = data.frame(hormon = c(0, 1)),
     col = c("steelblue", "firebrick"), main = "rpsurv: predicted survival")
```

`predict()` also returns `type = "hazard"`, `"cumhaz"`, or `"link"`
(the linear predictor), with `se.fit = TRUE` adding delta-method 95%
confidence limits for `"survival"`, `"cumhaz"`, and `"link"`.

## Time-varying effects (non-proportional hazards)

A covariate's *value* is fixed for a subject, but its *coefficient* can be
allowed to change over follow-up via its own spline in log time, requested
with `tve`:

```r
fit_tve <- rpsurv(Surv(rectime, censrec) ~ hormon, data = brcancer,
                   df = 4, tve = "hormon", tve.df = 2)
coef(fit_tve)
```

A likelihood-ratio test of `fit_tve` against `fit` (both maximum-likelihood
fits, nested models) tests the proportional-hazards assumption for
`hormon` directly.

This is distinct from a **time-varying covariate**, where the covariate's
*value itself* changes during follow-up. That is a data-representation
question, not a modeling one: split each subject's follow-up into
intervals over which the covariate is constant and supply counting-process
data via `Surv(start, stop, status)`. `rpsurv` detects the 3-column `Surv`
object and fits the left-truncated likelihood automatically, correctly
accounting for each interval's dependence on survival to its start time
under the covariate value of the *previous* interval. The same mechanism
also handles ordinary left truncation (delayed entry) when covariate
values don't change. See `vignette("rpsurv")` for a full worked example
and the underlying likelihood.

## Diagnostics

```r
residuals(fit, type = "coxsnell")   # also "martingale", "deviance"
coxsnell_plot(fit)                  # Cox-Snell residual plot
km_compare_plot(fit)                # fitted vs. Kaplan-Meier overlay
```

## Validation

Coefficients, standard errors, and log-likelihoods agree with
`rstpm2::stpm2` to 4-5 decimal places across all three scales (`hazard`,
`odds`, `normal`) and with `flexsurv::flexsurvspline` on the hazard scale.
See `tests/testthat/test-rpsurv.R` for the full parity suite, which also
covers left truncation.

## Benchmark vs. `rstpm2` and `flexsurv`

Wall-clock fit time on simulated Weibull-hazard data with two covariates,
reproducible via `data-raw/benchmark.R` (shipped with the package source,
not installed). `flexsurv::flexsurvspline` is materially slower per fit,
so it was only run up to n=2e4 to keep the benchmark tractable;
`rstpm2::stpm2` was run up to n=2e5.

| n | rpsurv (s) | rstpm2 (s) | flexsurvspline (s) |
|---:|---:|---:|---:|
| 1,000 | 0.031 | 0.017 | 0.121 |
| 5,000 | 0.037 | 0.058 | 0.636 |
| 20,000 | 0.080 | 0.238 | 2.330 |
| 50,000 | 0.165 | 1.290 | - |
| 100,000 | 0.247 | 2.410 | - |
| 200,000 | 0.413 | 4.912 | - |
| 500,000 | 1.182 | - | - |

At n=200,000, `rpsurv` is about 12x faster than `rstpm2::stpm2` and
continues to scale to n=500,000 (`rstpm2` was not run at that size). At
n=20,000, `rpsurv` is about 29x faster than `flexsurv::flexsurvspline`.
The advantage widens with n in both comparisons: `rpsurv` does O(n) work
per optimizer step in a single fused parallel pass over the log-likelihood
and its analytic gradient, computed together (`RcppParallel::parallelReduce`),
against repeated R-level likelihood evaluation (numerical gradients, or
gradients built from R-level matrix algebra) in the other two packages.
See `vignette("rpsurv")` for the full benchmark code and what else is
delegated to C++ (the spline basis itself, and closed-form starting
values from a bounded Kaplan-Meier subsample).

## Notes

- `scale` selects the transformation of the survival function the spline
  models: `"hazard"` (proportional hazards, the Royston-Parmar default),
  `"odds"` (proportional odds), or `"normal"` (probit).
- Boundary knots default to the min/max of log time among events, and
  interior knots to equally spaced centiles of log time among events
  (`default_knots()`), matching `rstpm2`/`flexsurv` conventions.
- `SystemRequirements: GNU make` (for `RcppParallel`); no other external
  dependencies.
