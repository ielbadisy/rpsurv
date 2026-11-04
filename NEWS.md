# rpsurv 0.8.2

* `summary()` no longer shows `exp(coef)` for the baseline spline terms: it has no interpretation, and a tiny `exp(intercept)` pushed the whole table into scientific notation. Tables follow the session's `digits` and `scipen` options; log-likelihood, AIC and BIC are printed with two decimals in `print()` and `summary()`.

* `predict()`, `rmst_diff()`, `plot()` and the other functions taking
  `newdata` / `newdata0` now accept covariates on the original variables
  of the formula, like `coxph()` and `stpm2()`: factors or character
  values by their levels, untransformed numeric variables for terms such
  as `log(karno)` or `poly(age, 2)`, and interactions. Previously a factor
  had to be given as its treatment-contrast indicator columns (for example
  `celltypeadeno`), otherwise `predict()` stopped with "`newdata` is
  missing covariate(s)". The fit now stores its terms, factor levels and
  contrasts; levels not seen when fitting give a clear error, and data
  holding the indicator columns is still accepted. Naming a factor in
  `tve` or `nle` now lists its indicator columns in the error.
* `km_compare_plot()` compared Kaplan-Meier with the survival curve at
  the mean covariate profile, a different estimand that is biased when
  covariates vary. The model curve is now the population-averaged
  survival of the group (mean of the individual predicted curves), and
  for counting-process data (left truncation, time-varying covariates)
  the cumulated average hazard of the rows at risk, which is what the
  left-truncated Kaplan-Meier estimates. `by` failed for factors ("group
  length is 0"); it now takes any variable of the formula or a grouping
  vector (for example `cut(age, 3)`), gives each group its own colour and
  a legend saying which group is which, and the plotted coordinates are
  returned invisibly.
* `plot.rpsurv()`: the y axis now covers the confidence band (it was set
  from the point estimates only, so the band was clipped); survival is
  drawn on [0, 1]. The default profile no longer averages factor
  indicator columns into a non-existent "mean" level: numeric variables
  are at their mean and factors at their most common level. Several
  curves get a legend labelled from the `newdata` rows, placed in the
  emptiest corner (`legend` argument). `type = "hazard"` now also draws
  its confidence band.
* `plot.tvecurve()`: the default y label said "per unit" of the
  covariate for `type = "hr"`, but the hazard ratio is for `ref + 1`
  versus `ref` and depends on `ref` under a time-varying effect; the label
  now says so. For `type = "link"` the label names the model scale (log
  cumulative hazard ratio, log cumulative odds ratio, or probit-scale
  effect) instead of a generic "beta(t)".
* `coxsnell_plot()` now uses the Nelson-Aalen estimator as documented
  (it used `-log` of Kaplan-Meier, infinite when the last residual is an
  event), and draws both axes on the same range from 0 so the y = x line
  is the diagonal.
* `plot.rpsurv()`, `km_compare_plot()` and `coxsnell_plot()` accept `xlab`
  and `ylab` in `...` (previously an error, "formal argument matched by
  multiple actual arguments").

* `summary()` and `confint()` give NA standard errors with one clear warning when the vcov is not positive definite, instead of R's generic "NaNs produced".

* Fixed non-positive-definite variance-covariance matrices (NaN standard
  errors, "NaNs produced" in `summary()`) for fits with a time-varying
  effect of a covariate on a large scale, e.g. `tve = "karno"` on
  `veteran`. The covariate times the cubic spline in log time made the
  optimisation badly conditioned, so the optimiser stopped early and the
  finite-difference Hessian was indefinite. `rpsurv()` now optimises over
  an orthogonalised linear reparameterisation of the design (from its QR
  decomposition) and maps the coefficients and variance-covariance matrix
  back exactly. Standard errors of `tve` fits on continuous covariates now
  agree with `rstpm2::stpm2()`; earlier finite standard errors for such
  fits could also be inaccurate. Point estimates are unchanged up to
  optimiser tolerance. Settings in `control` such as `parscale` now refer
  to the internal parameterisation.
* Factor covariates are now coded by treatment contrasts (one column per
  non-reference level). Previously the first factor in the formula got an
  indicator column for every level, which is collinear with the model
  intercept: the model was not identified and its standard errors were
  NaN, without a warning. Factor fits now match `rstpm2::stpm2()`.
* `predict(type = "cumhaz", se.fit = TRUE)` returned the confidence
  limits swapped (`lower` above `upper`). The cumulative hazard and its
  limits are now computed directly from the linear predictor, so they no
  longer overflow to `Inf` when the survival probability underflows to 0
  (also used by `residuals()`).
* `rpsurv()` no longer fails silently. It now warns when `optim()` does
  not converge (`convergence != 0`), when the optimiser stops short of a
  stationary point, when the design matrix is rank deficient (naming the
  aliased columns), when the Hessian is singular or not finite (standard
  errors `NA`), when the Hessian is not positive definite (one warning; the
  variance-covariance matrix is then not positive definite and
  `summary()` / `confint()` report the affected standard errors as `NA`), and when the final estimates sit on the
  likelihood penalty for a decreasing cumulative hazard. Tied centiles
  that merge spline knots now warn, like knots dropped at a boundary.
* `rpsurv()` now stops with an informative error for data without events,
  non-integer or non-positive `df`, `tve.df` or `nle.df` (previously
  truncated silently), and duplicated user `knots`.
* `print()` and `summary()` report observations dropped because of missing
  values, as `print.coxph()` does.
* `predict()`, `tvecurve()` and `rmst_diff()` now reject non-positive or
  non-finite `times` / `tau` (time 0 previously gave `NaN` silently).
  `predict()`, `tvecurve()`, `nlecurve()` and `rmst_diff()` warn once,
  with the reason, when the fitted hazard is not positive at some
  evaluation points (hazard-based estimates are then negative or `NaN`),
  when the variance-covariance matrix is unavailable, or when standard
  errors are `NaN` because it is not positive definite. `residuals()`
  warns about negative or non-finite Cox-Snell residuals.
* New tests validate `nle` fits, alone and combined with `tve`, on the
  hazard, odds and probit scales against `rstpm2::stpm2()` with the same
  natural-spline model space (log-likelihood, coefficients and full
  variance-covariance matrix), and check `nlecurve()` and `tvecurve()`
  standard errors against a numerical delta method.

# rpsurv 0.8.1

* `tvecurve()` and `nlecurve()` now return data frames of class
  `"tvecurve"` and `"nlecurve"`, which keep the covariate, reference value
  and scale as attributes.
* New `plot()` methods for both: the hazard ratio (or log effect) against
  time or against the covariate, with its 95% confidence band, a reference
  line at no effect and a log y axis by default. `add = TRUE` overlays
  curves.

# rpsurv 0.8.0

* New `nle` and `nle.df` arguments to `rpsurv()`: a smooth non-linear
  effect g(x) for numeric covariates, modelled as a restricted cubic
  spline in the covariate (knots at its range and equally spaced
  centiles, basis rescaled to [0, 1] for conditioning). It is specified
  separately from the time-varying effect (`tve`, `tve.df`), and a
  covariate may have both, giving g(x) + beta(t) x. `nle.df = 1`
  reproduces the linear fit exactly.
* `tve.df` and `nle.df` accept a single value or a vector named by
  covariate.
* New `nlecurve()`: g(x) - g(ref) with a delta-method 95% confidence
  interval (the log hazard ratio versus `ref` on the hazard scale when
  the covariate is not also time-varying).
* New `tvecurve()`: the time-varying effect, either as the log
  instantaneous hazard ratio per unit increase (`type = "hr"`, matching
  `predict(type = "hr")`) or as the coefficient function beta(t) on the
  linear-predictor scale (`type = "link"`). For a covariate in both `nle`
  and `tve`, `tvecurve()` and `nlecurve()` return the two components
  separately, and their sum equals the total log hazard ratio.

# rpsurv 0.7.1

* CRAN resubmission: added a `<doi:10.1002/sim.1203>` link for the
  Royston and Parmar (2002) reference in the `DESCRIPTION` Description
  field, and added `\value` documentation for `coxsnell_plot()`,
  `km_compare_plot()`, and `plot.rpsurv()`.

# rpsurv 0.7.0

* Added `predict.rpsurv(type = "hr")`: the model-implied instantaneous
  hazard ratio between a contrast covariate profile (`newdata`) and a
  reference profile (`newdata0`), computed as the ratio of two
  `type = "hazard"` predictions. This is the correct contrast under a
  time-varying effect (`tve`), where `exp(eta1 - eta0)` is not the
  instantaneous hazard ratio in general. Matches `exp(coef)` exactly
  under proportional hazards and agrees with `rstpm2::predict(...,
  type = "hr")` under `tve` to within Monte Carlo/spline-basis tolerance.
  `se.fit` is not yet supported for `type = "hr"` (nor for `"hazard"`).

# rpsurv 0.6.1

## CRAN submission prep

* Bundled three datasets, each chosen for a specific package feature
  rather than as a generic example: `brcancer` (baseline fitting example,
  identical to `rstpm2::brcancer`), `veteran` (a well-known
  non-proportional-hazards effect, illustrating `tve`, identical to
  `survival::veteran`), and `heart` (Stanford heart transplant data,
  already in counting-process format with a time-varying covariate,
  identical to `survival::heart`). Removes the need for `rstpm2` to be
  installed just to run the README/vignette examples.
* Added a detailed `README.md` with real executed output (fit/summary/
  predict, time-varying effects, diagnostics, and the `rstpm2`/`flexsurv`
  speed benchmark table).
* Added a GitHub Actions `R-CMD-check.yaml` workflow and matching README
  badges.
* Fixed `.Rbuildignore`: `data-raw/` (containing a cached benchmark `.rds`)
  was shipping in the source tarball as a non-standard top-level
  directory; also ignored `cran-comments.md` and `*.Rcheck`.
* Added `cran-comments.md` for the first CRAN submission.

# rpsurv 0.6.0

* Standardized author name casing to "El Badisy" (title case) in
  `Authors@R` and citation-facing fields, matching the rest of the
  package suite.

# rpsurv 0.5.0

## Performance

* Delegated the restricted cubic spline basis (`rcs_basis()`) to a C++
  implementation (`rcs_basis_cpp()`), rather than R's vectorised-but-
  interpreted arithmetic.
* Cached the fused log-likelihood/gradient C++ call keyed on the
  parameter vector, since `optim(method = "BFGS")` otherwise evaluates
  `fn` and `gr` separately at every trial point, doubling data passes for
  the same parameter vector.
* Closed-form starting values are now computed from a bounded (default
  20,000-row) subsample of the Kaplan-Meier curve, avoiding an
  O(n log n) sort of the full data purely for initialization.
* Added a reproducible speed benchmark (`data-raw/benchmark.R`) against
  `rstpm2::stpm2` and `flexsurv::flexsurvspline`, plus a full vignette
  covering the model, API, and benchmark.

# rpsurv 0.4.0

* Added support for left-truncated / counting-process data
  (`Surv(start, stop, status)`), enabling genuine time-varying covariates
  (as opposed to time-varying *effects*, see `tve`): each interval
  contributes `log S(stop) - log S(start)`, correctly conditioning on
  survival to `start` under the covariate values of the *previous*
  interval. The same mechanism handles ordinary left truncation (delayed
  entry) when covariate values don't change.
* Added deviance residuals and `km_compare_plot()` (fitted vs.
  Kaplan-Meier calibration diagnostic).
* Renamed `tvc` to `tve` (time-varying effect) throughout the fitting,
  predict, and summary API, to avoid confusion with the new
  time-varying-*covariate* support above (a different, data-representation-
  level concept).

# rpsurv 0.3.0

* Added `predict.rpsurv()` (`type = "survival"`/`"hazard"`/`"cumhaz"`/
  `"link"`, with optional delta-method confidence limits) and
  `plot.rpsurv()`.
* Added `print`/`summary`/`coef`/`vcov`/`logLik`/`AIC`/`BIC`/`confint`
  methods, matching `coxph()`'s output style: interpretable covariate
  effects (hazard/odds ratios, Wald tests) are reported separately from
  the baseline spline's nuisance coefficients.
* Added Cox-Snell and martingale residuals plus `coxsnell_plot()`.
* Validated against `rstpm2::stpm2` and `flexsurv::flexsurvspline` on
  `rstpm2::brcancer` across all three link scales (hazard, odds, normal);
  coefficients, standard errors, and log-likelihoods agree to 4-5 decimal
  places.

# rpsurv 0.2.0

* Added the parallel C++ log-likelihood and its analytic gradient
  (`RcppParallel::parallelReduce`), and `rpsurv()`, the main model-fitting
  function, with support for a time-varying effect (`tve`) via its own
  spline in log time multiplying the covariate.

# rpsurv 0.1.0

* Initial implementation: restricted cubic spline basis in log time
  (`rcs_basis()`) and the Royston-Parmar design matrix builder, following
  Royston and Parmar (2002).
