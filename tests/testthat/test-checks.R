test_that("non-convergence of optim() gives a warning", {
  expect_warning(
    rpsurv(survival::Surv(time, status) ~ karno, data = veteran, control = list(maxit = 3)),
    "did not converge"
  )
})

test_that("a rank-deficient design and a singular Hessian are reported", {
  d <- veteran
  d$k2 <- 2 * d$karno
  w <- character()
  fit <- withCallingHandlers(
    rpsurv(survival::Surv(time, status) ~ karno + k2, data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("rank deficient.*k2", w)))
  expect_true(any(grepl("singular", w)))
  expect_true(all(is.na(vcov(fit))))
})

test_that("a non-positive-definite Hessian warns once", {
  d <- veteran
  d$start <- pmin(d$time / 2, 5)
  w <- character()
  fit <- withCallingHandlers(
    rpsurv(survival::Surv(start, time, status) ~ karno, data = d, tve = "karno"),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("not positive definite", w), 1L)
})

test_that("well-behaved fits give no warnings", {
  expect_silent(rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran, tve = "karno"))
  expect_silent(rpsurv(survival::Surv(time, status) ~ trt + karno + age, data = veteran,
                       nle = "age", scale = "odds"))
  for (k in 1:6) {
    expect_silent(rpsurv(survival::Surv(time, status) ~ karno, data = veteran, df = k))
  }
})

test_that("invalid inputs give informative errors", {
  d <- veteran
  d$status <- 0
  expect_error(rpsurv(survival::Surv(time, status) ~ karno, data = d), "no events")
  expect_error(rpsurv(survival::Surv(time, status) ~ karno, data = veteran, df = 2.5), "whole number")
  expect_error(rpsurv(survival::Surv(time, status) ~ karno, data = veteran, df = 0), "whole number")
  expect_error(rpsurv(survival::Surv(time, status) ~ karno, data = veteran, tve = "karno",
                      tve.df = 1.5), "whole number")
  expect_error(rpsurv(survival::Surv(time, status) ~ karno, data = veteran,
                      knots = log(c(1, 50, 50, 999))), "duplicated")
})

test_that("tied centiles that merge knots are reported", {
  expect_warning(k <- default_knots(c(rep(0, 10), rep(1, 30), rep(2, 30), 3), 4), "tied")
  expect_equal(k, c(0, 1, 2, 3))
})

test_that("print and summary report observations dropped for missingness", {
  d <- veteran
  d$karno[1:5] <- NA
  fit <- rpsurv(survival::Surv(time, status) ~ karno, data = d)
  expect_equal(fit$n, nrow(d) - 5L)
  expect_output(print(fit), "5 observations deleted due to missingness")
  expect_output(print(summary(fit)), "5 observations deleted due to missingness")
})
