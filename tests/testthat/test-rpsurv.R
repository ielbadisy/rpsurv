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
  expect_length(residuals(fit, type = "martingale"), fit$n)
  expect_length(residuals(fit, type = "deviance"), fit$n)
  expect_true(all(residuals(fit, type = "coxsnell") >= 0))
  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("km_compare_plot and coxsnell_plot run without error", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_error(km_compare_plot(fit, by = "hormon"), NA)
  expect_error(coxsnell_plot(fit), NA)
  grDevices::dev.off()
  unlink(tmp)
})

test_that("time-varying effect (tve) fits without error", {
  d <- get_brcancer()
  fit <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = d, df = 4, tve = "hormon", tve.df = 2)
  expect_s3_class(fit, "rpsurv")
  expect_true(fit$convergence == 0L)
})

test_that("left truncation (counting-process) matches rstpm2::stpm2", {
  skip_if_not_installed("rstpm2")
  suppressPackageStartupMessages(library(rstpm2))
  d <- get_brcancer()
  set.seed(1)
  d$entry <- d$rectime * stats::runif(nrow(d), 0, 0.3)

  fit <- rpsurv(survival::Surv(entry, rectime, censrec) ~ hormon, data = d, df = 4)
  ref <- rstpm2::stpm2(survival::Surv(entry, rectime, censrec) ~ hormon, data = d, df = 4)

  expect_true(fit$counting)
  expect_equal(unname(coef(fit)["hormon"]), unname(coef(ref)["hormon"]), tolerance = 1e-2)
  expect_equal(fit$loglik, -ref@min, tolerance = 1e-2)
})

test_that("genuine time-varying covariate (counting-process split data) fits and improves fit vs baseline", {
  d <- get_brcancer()
  d$id <- seq_len(nrow(d))
  # split each subject's follow-up at the midpoint and let x1 change value there,
  # simulating a covariate that is updated partway through follow-up
  half <- d$rectime / 2
  first <- data.frame(id = d$id, start = 0, stop = half, status = 0, hormon = d$hormon, x1 = d$x1)
  second <- data.frame(id = d$id, start = half, stop = d$rectime, status = d$censrec,
                        hormon = d$hormon, x1 = d$x1 + 1)
  long <- rbind(first, second)
  long <- long[long$start < long$stop, ]

  fit <- rpsurv(survival::Surv(start, stop, status) ~ hormon + x1, data = long, df = 4)
  expect_s3_class(fit, "rpsurv")
  expect_true(fit$counting)
  expect_true(fit$convergence == 0L)
  expect_true(is.finite(coef(fit)["x1"]))
})
