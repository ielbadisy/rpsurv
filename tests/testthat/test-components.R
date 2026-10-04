sim_nle <- function(n = 1500, seed = 11) {
  set.seed(seed)
  x1 <- rbinom(n, 1, 0.5)
  x2 <- runif(n, 0, 3)
  g <- 0.5 * (x2 - 1.5)^2
  t_ev <- rweibull(n, shape = 1.2, scale = 1 / (0.3 * exp(0.7 * x1 + g))^(1 / 1.2))
  t_c <- rexp(n, 0.1)
  data.frame(time = pmin(t_ev, t_c), status = as.integer(t_ev <= t_c), x1 = x1, x2 = x2)
}

test_that("nle.df = 1 reproduces the linear fit", {
  d <- sim_nle()
  f_lin <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d)
  f_nle <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, nle = "x2", nle.df = 1)
  expect_equal(f_nle$loglik, f_lin$loglik, tolerance = 1e-6)
  expect_equal(nlecurve(f_nle, "x2", at = 2, ref = 1)$est,
               unname(coef(f_lin)["x2"]), tolerance = 1e-4)
})

test_that("nlecurve recovers a known non-linear effect and matches predict(type = 'hr')", {
  d <- sim_nle(n = 4000)
  f <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, nle = "x2", nle.df = 4)
  at <- c(0.5, 1, 2, 2.5)
  g <- nlecurve(f, "x2", at = at, ref = 1.5)
  expect_equal(g$est, 0.5 * (at - 1.5)^2, tolerance = 0.25)
  expect_true(all(g$lower < g$est & g$est < g$upper))
  hr <- predict(f, newdata = data.frame(x1 = 0, x2 = at), newdata0 = data.frame(x1 = 0, x2 = rep(1.5, 4)),
                times = 1, type = "hr")
  expect_equal(log(hr$est), g$est, tolerance = 1e-8)
})

test_that("tvecurve matches predict(type = 'hr') and the link-scale coefficient", {
  d <- sim_nle()
  f <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, tve = "x1", nle = "x2")
  tt <- c(0.5, 1, 2)
  tv <- tvecurve(f, "x1", times = tt)
  hr <- predict(f, newdata = data.frame(x1 = 1, x2 = 0), newdata0 = data.frame(x1 = 0, x2 = 0),
                times = tt, type = "hr", se.fit = TRUE)
  expect_equal(tv$est, log(hr$est), tolerance = 1e-8)
  expect_equal(tv$upper - tv$est, log(hr$upper) - log(hr$est), tolerance = 1e-8)
  ln <- tvecurve(f, "x1", times = tt, type = "link")
  b <- coef(f)
  manual <- b["x1"] + rpsurv:::rcs_basis(log(tt), f$tve_knots$x1) %*% b[grep("^x1:s\\(logt\\)", names(b))]
  expect_equal(ln$est, as.numeric(manual), tolerance = 1e-10)
})

test_that("a covariate in both nle and tve decomposes into the two components", {
  d <- sim_nle()
  f <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, tve = "x2", nle = "x2")
  tt <- c(0.5, 1, 2)
  total <- predict(f, newdata = data.frame(x1 = 0, x2 = 2), newdata0 = data.frame(x1 = 0, x2 = 1),
                   times = tt, type = "hr")
  parts <- tvecurve(f, "x2", times = tt, ref = 1)$est + nlecurve(f, "x2", at = 2, ref = 1)$est
  expect_equal(log(total$est), parts, tolerance = 1e-8)
})

test_that("per-covariate df and input checks", {
  d <- sim_nle()
  f <- rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, tve = "x1", tve.df = c(x1 = 2),
              nle = "x2", nle.df = c(x2 = 3))
  expect_length(grep("^x1:s\\(logt\\)", names(coef(f))), 2)
  expect_length(grep("^x2:s\\(x\\)", names(coef(f))), 3)
  expect_error(rpsurv(survival::Surv(time, status) ~ x1, data = d, nle = "x2"), "not found")
  expect_error(rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, nle = "x2", nle.df = c(x3 = 2)), "named")
  expect_error(rpsurv(survival::Surv(time, status) ~ x1 + x2, data = d, nle = "x1"), "distinct")
  expect_error(nlecurve(f, "x1"), "nle")
  expect_error(tvecurve(f, "x2"), "tve")
})

test_that("nle predictions work for left-truncated data and new covariate values", {
  d <- sim_nle()
  d$start <- pmin(d$time / 3, 0.2)
  f <- rpsurv(survival::Surv(start, time, status) ~ x1 + x2, data = d, nle = "x2")
  p <- predict(f, newdata = data.frame(x1 = 0, x2 = c(-0.5, 3.5)), times = 1)
  expect_true(all(is.finite(p$est)) && all(p$est > 0 & p$est < 1))
})
