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
