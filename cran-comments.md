# CRAN submission comments: rpsurv 0.7.1

## Resubmission

This is a resubmission. In response to the CRAN reviewer's comments on
0.6.1:

* Added a reference link to the Description field: the Royston and Parmar
  (2002) method reference is now given as `<doi:10.1002/sim.1203>`, with
  no space after `doi:`.
* Added `\value` sections to the exported `coxsnell_plot()`,
  `km_compare_plot()` and `plot.rpsurv()` help pages, documenting the
  class and meaning of the returned object (or that the function is
  called for its plotting side effect).

## Test environments



* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 ERRORs | 0 WARNINGs | 2-3 NOTEs

* "New submission", standard for a first CRAN submission.
* "unable to verify current time" appears intermittently, a local
  clock-check note unrelated to the package.
* "Compilation used the following non-portable flag(s): -mno-omit-leaf-frame-pointer".
  This comes from the local Ubuntu-packaged R toolchain's own `CFLAGS`
  (confirmed via `R CMD config CFLAGS`), not from the package's own
  `src/Makevars` (which already filters this flag out of `CXX14FLAGS`).
  Not expected to reproduce on CRAN's own build machines.

`R CMD build --compact-vignettes=both` was used to compact the vignette PDF
(447Kb to 175Kb), which otherwise triggers a "checking sizes of PDF files
under 'inst/doc'" WARNING under `--as-cran`.

## Downstream dependencies

There are currently no downstream dependencies for this package, as it has
not previously been published on CRAN.
