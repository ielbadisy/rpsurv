# CRAN submission comments: rpsurv 0.6.1

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 ERRORs | 0 WARNINGs | 3 NOTEs

* "New submission", standard for a first CRAN submission.
* "unable to verify current time", a local clock-check note unrelated to
  the package.
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
