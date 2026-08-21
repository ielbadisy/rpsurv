## Bundles survival::heart (Stanford heart transplant study, Crowley &
## Hu 1977), already in counting-process (start, stop, event) format with
## a time-varying `transplant` covariate. Used to illustrate rpsurv's
## left-truncated / counting-process likelihood for genuine time-varying
## covariates (Surv(start, stop, status)), without needing to construct a
## synthetic example.
heart <- survival::heart
usethis::use_data(heart, overwrite = TRUE)
