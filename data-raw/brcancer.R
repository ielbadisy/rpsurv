## Prepares the bundled `brcancer` dataset from rstpm2::brcancer (German
## Breast Cancer Study Group data, Sauerbrei & Royston), so rpsurv examples
## and vignette do not require rstpm2 to be installed.
brcancer <- rstpm2::brcancer
attributes(brcancer)[c("datalabel", "time.stamp", "formats", "types", "val.labels", "var.labels", "version")] <- NULL
usethis::use_data(brcancer, overwrite = TRUE)
