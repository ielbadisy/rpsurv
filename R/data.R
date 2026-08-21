#' German breast cancer data
#'
#' Recurrence-free survival for 686 women in a randomized trial of
#' hormonal therapy for breast cancer (Sauerbrei and Royston), used
#' throughout the package's documentation and vignette as the baseline
#' `rpsurv()` fitting example. Identical to `rstpm2::brcancer`, bundled
#' here so examples do not require `rstpm2` to be installed.
#'
#' @format A data frame with 686 rows and 14 variables:
#' \describe{
#'   \item{id}{subject identifier}
#'   \item{hormon}{hormonal therapy (0 = no, 1 = yes)}
#'   \item{x1}{age, years}
#'   \item{x2}{menopausal status}
#'   \item{x3}{tumour size, mm}
#'   \item{x4}{tumour grade}
#'   \item{x5}{number of positive lymph nodes}
#'   \item{x6}{progesterone receptor, fmol}
#'   \item{x7}{estrogen receptor, fmol}
#'   \item{rectime}{recurrence-free survival time, days}
#'   \item{censrec}{event indicator (1 = recurrence/death, 0 = censored)}
#'   \item{x4a}{tumour grade >= 2}
#'   \item{x4b}{tumour grade == 3}
#'   \item{x5e}{exp(-0.12 * x5)}
#' }
#' @source \url{https://www.stata-press.com/data/r11/brcancer.dta}, as
#'   redistributed by the \pkg{rstpm2} package.
"brcancer"

#' Veterans' Administration lung cancer trial
#'
#' A randomized trial comparing two treatment regimens for lung cancer
#' (Kalbfleisch and Prentice). `celltype` is a well-known example of a
#' non-proportional-hazards effect: its hazard ratio changes materially
#' over follow-up, which makes this dataset a natural illustration of
#' `rpsurv()`'s `tve` (time-varying effect) argument. Identical to
#' `survival::veteran`, bundled here under the same name for convenience.
#'
#' @format A data frame with 137 rows and 8 variables:
#' \describe{
#'   \item{trt}{treatment (1 = standard, 2 = test)}
#'   \item{celltype}{tumour cell type: squamous, smallcell, adeno, large}
#'   \item{time}{survival time, days}
#'   \item{status}{event indicator (1 = dead, 0 = censored)}
#'   \item{karno}{Karnofsky performance score (0-100)}
#'   \item{diagtime}{months from diagnosis to randomization}
#'   \item{age}{age, years}
#'   \item{prior}{prior therapy (0 = no, 10 = yes)}
#' }
#' @source Kalbfleisch, J. and Prentice, R. (1980) *The Statistical
#'   Analysis of Failure Time Data*. Wiley, New York, as redistributed by
#'   the \pkg{survival} package.
"veteran"

#' Stanford heart transplant data
#'
#' Survival of patients on the waiting list for the Stanford heart
#' transplant program (Crowley and Hu, 1977), already in counting-process
#' (`start`, `stop`, `event`) format with a time-varying `transplant`
#' covariate: a patient's row is split at the moment of transplantation,
#' if it occurred. This makes it a natural, real (rather than
#' constructed) illustration of `rpsurv()`'s support for genuine
#' time-varying covariates via `Surv(start, stop, status)`. Identical to
#' `survival::heart`, bundled here under the same name for convenience.
#'
#' @format A data frame with 172 rows and 8 variables:
#' \describe{
#'   \item{start}{interval start time, days}
#'   \item{stop}{interval stop time, days}
#'   \item{event}{event indicator at `stop` (1 = death, 0 = censored)}
#'   \item{age}{age minus 48 years}
#'   \item{year}{year of acceptance into the program, in years since 1967}
#'   \item{surgery}{prior bypass surgery (0 = no, 1 = yes)}
#'   \item{transplant}{transplant status during this interval (0 = no, 1 = yes)}
#'   \item{id}{patient identifier}
#' }
#' @source Crowley, J. and Hu, M. (1977) Covariance analysis of heart
#'   transplant survival data. *Journal of the American Statistical
#'   Association*, 72, 27-36, as redistributed by the \pkg{survival}
#'   package.
"heart"
