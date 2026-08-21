## Bundles survival::veteran (Veterans' Administration lung cancer trial,
## Kalbfleisch & Prentice). Classic non-proportional-hazards example: the
## `celltype` effect is well known not to satisfy proportional hazards,
## which makes it a natural illustration of rpsurv's `tve` argument.
veteran <- survival::veteran
usethis::use_data(veteran, overwrite = TRUE)
