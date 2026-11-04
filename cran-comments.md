## Update submission (0.8.1 -> 0.8.2)

This is a bug-fix release soon after 0.8.1. It fixes wrong results that
0.8.1 returns without any warning:

* Factor covariates: the first factor in the formula got an indicator for
  every level, collinear with the intercept, so the model was not
  identified and the standard errors were NaN. Factors now use treatment
  contrasts.
* Time-varying effects of continuous covariates (`tve`): the optimisation
  was badly conditioned, so the variance-covariance matrix could be
  indefinite (NaN standard errors) or inaccurate. The model is now fitted
  on an orthogonalised reparameterisation and mapped back exactly.
* `predict(type = "cumhaz", se.fit = TRUE)` returned swapped confidence
  limits.
* `predict()` and the plot functions now accept factor covariates by
  their levels in `newdata`; `km_compare_plot()` compared Kaplan-Meier
  with the curve at the mean covariate profile instead of the average
  predicted survival, and its `by` argument failed for factors.

Fits are now validated against `rstpm2::stpm2()` in the tests
(coefficients, full variance-covariance matrix, log-likelihood). Fitting
and prediction problems that were silent before now give warnings or
errors. See NEWS.md for details.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

* "Days since last update: 1": a bug-fix release, see above.

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This comes from the local Ubuntu R
  toolchain (`R CMD config CFLAGS`), not from the package's `src/Makevars`.

The tarball was built with `R CMD build --compact-vignettes=both` to
compact the vignette PDF.

## Reverse dependencies

None on CRAN.
