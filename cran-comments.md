## Update submission (0.7.1 -> 0.8.1)

See NEWS.md for details.

* New `nle` and `nle.df` arguments to `rpsurv()`: a smooth non-linear
  covariate effect (restricted cubic spline), specified separately from
  the time-varying effect (`tve`). `nle.df = 1` reproduces the linear fit
  exactly.
* `tve.df` and `nle.df` accept a single value or a named vector.
* New `nlecurve()` and `tvecurve()` with delta-method 95% confidence
  intervals, and `plot()` methods for both.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This comes from the local Ubuntu R
  toolchain (`R CMD config CFLAGS`), not from the package's `src/Makevars`.

The tarball was built with `R CMD build --compact-vignettes=both` to
compact the vignette PDF.

## Reverse dependencies

None on CRAN.
