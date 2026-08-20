skip_rstpm2 <- !requireNamespace("rstpm2", quietly = TRUE)

get_brcancer <- function() {
  d <- get(data("brcancer", package = "rstpm2", envir = environment()))
  d$hormon <- as.numeric(d$hormon)
  d
}

test_that("basic PH fit recovers a sensible model on brcancer", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, scale = "hazard")
  expect_s3_class(fit, "rpsurv")
  expect_equal(fit$convergence, 0L)
  expect_true(fit$loglik > -3000)
  expect_true(is.finite(sqrt(diag(fit$vcov))["hormon"]))
})

test_that("PH scale matches rstpm2::stpm2 coefficients and log-likelihood", {
  skip_if_not_installed("rstpm2")
  suppressPackageStartupMessages(library(rstpm2)) # stpm2() calls internal gsm() unqualified
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, scale = "hazard")
  ref <- rstpm2::stpm2(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4)

  expect_equal(unname(coef(fit)["hormon"]), unname(coef(ref)["hormon"]), tolerance = 1e-3)
  expect_equal(fit$loglik, -ref@min, tolerance = 1e-4)
})

test_that("PO and probit scales match rstpm2::stpm2", {
  skip_if_not_installed("rstpm2")
  suppressPackageStartupMessages(library(rstpm2)) # stpm2() calls internal gsm() unqualified
  d <- get_brcancer()

  fit_po <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, scale = "odds")
  ref_po <- rstpm2::stpm2(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, link.type = "PO")
  expect_equal(unname(coef(fit_po)["hormon"]), unname(coef(ref_po)["hormon"]), tolerance = 1e-3)

  fit_pr <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, scale = "normal")
  ref_pr <- rstpm2::stpm2(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, link.type = "probit")
  expect_equal(unname(coef(fit_pr)["hormon"]), unname(coef(ref_pr)["hormon"]), tolerance = 1e-3)
})

test_that("predict returns monotone decreasing survival curves", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4)
  pred <- predict(fit, newdata = data.frame(hormon = 0), type = "survival")
  expect_true(all(diff(pred$est) <= 1e-8))
  expect_true(all(pred$est >= 0 & pred$est <= 1))
})

test_that("residuals and print/summary methods run without error", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4)
  expect_length(residuals(fit, type = "coxsnell"), fit$n)
  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("time-varying effect (tvc) fits without error", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, tvc = "hormon", tvc.df = 2)
  expect_s3_class(fit, "rpsurv")
  expect_true(fit$convergence == 0L)
})
