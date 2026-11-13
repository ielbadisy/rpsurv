test_that("newdata takes factor and character covariates on their original variables", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype + karno, data = veteran)
  nd <- data.frame(trt = c(1, 2), celltype = c("adeno", "squamous"), karno = c(60, 70))
  p <- predict(fit, newdata = nd, times = c(30, 100), se.fit = TRUE)
  expect_equal(nrow(p), 4L)
  expect_true(all(is.finite(p$est)))

  dummies <- data.frame(trt = c(1, 2), celltypesmallcell = 0, celltypeadeno = c(1, 0),
                        celltypelarge = 0, karno = c(60, 70))
  expect_equal(predict(fit, newdata = dummies, times = c(30, 100), se.fit = TRUE), p)
  nd_f <- nd
  nd_f$celltype <- factor(nd_f$celltype, levels = levels(veteran$celltype))
  expect_equal(predict(fit, newdata = nd_f, times = c(30, 100), se.fit = TRUE), p)

  expect_error(predict(fit, newdata = data.frame(trt = 1, celltype = "foo", karno = 60), times = 10),
               "level\\(s\\) of celltype not seen")
  expect_error(predict(fit, newdata = data.frame(trt = 1, karno = 60), times = 10),
               "missing covariate\\(s\\): celltype")

  hr <- predict(fit, newdata = nd[1, ], newdata0 = transform(nd[1, ], celltype = "squamous"),
                times = 50, type = "hr")
  expect_equal(hr$est, unname(exp(coef(fit)["celltypeadeno"])), tolerance = 1e-8)
  rd <- rmstdiff(fit, nd[1, ], transform(nd[1, ], celltype = "squamous"), tau = 200)
  expect_true(rd$est < 0)
})

test_that("factor and transformed-term predictions match rstpm2", {
  # reference: rstpm2 1.x, stpm2(Surv(time, status) ~ celltype + log(karno),
  # data = veteran, df = 4), predict(type = "surv", se.fit = TRUE) at t = 100
  fit <- rpsurv(survival::Surv(time, status) ~ celltype + log(karno), data = veteran)
  nd <- data.frame(celltype = c("adeno", "large"), karno = c(40, 80))
  p <- predict(fit, newdata = nd, times = 100, se.fit = TRUE)
  expect_equal(p$est, c(0.10794339409, 0.63369200406), tolerance = 1e-5)
  expect_equal(p$lower, c(0.03683529147, 0.49769617333), tolerance = 1e-4)
  expect_equal(p$upper, c(0.22287428104, 0.74211280388), tolerance = 1e-4)
})

test_that("interactions with factors are coded from the original variables", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt * celltype, data = veteran)
  p <- predict(fit, newdata = data.frame(trt = 2, celltype = "large"), times = 50, type = "link")
  b <- coef(fit)
  lp <- predict(fit, newdata = data.frame(trt = 0, celltype = "squamous"), times = 50, type = "link")$est +
    2 * b[["trt"]] + b[["celltypelarge"]] + 2 * b[["trt:celltypelarge"]]
  expect_equal(p$est, lp, tolerance = 1e-10)
})

test_that("tve on a factor name points to its indicator columns", {
  expect_error(rpsurv(survival::Surv(time, status) ~ celltype, data = veteran, tve = "celltype"),
               "celltypeadeno")
})

test_that("plot.rpsurv returns the plotted predictions, band included in the y range", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype, data = veteran)
  nd <- data.frame(trt = 1, celltype = c("squamous", "adeno"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (type in c("survival", "cumhaz", "hazard")) {
    p <- plot(fit, newdata = nd, type = type)
    expect_equal(p, predict(fit, newdata = nd, type = type, se.fit = TRUE))
    expect_false(is.null(p$lower))
    yr <- graphics::par("usr")[3:4]
    expect_true(yr[1] <= min(p$lower) && yr[2] >= max(p$upper))
  }
  expect_silent(plot(fit, xlab = "Days", main = "x", legend = FALSE))
  expect_equal(rpsurv:::rp_profile_labels(nd), c("celltype = squamous", "celltype = adeno"))
})

test_that("plot.rpsurv default profile uses the most common level of a factor", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype, data = veteran)
  prof <- rpsurv:::rp_default_profile(fit)
  expect_equal(prof$trt, mean(veteran$trt))
  expect_equal(as.character(prof$celltype), names(which.max(table(veteran$celltype))))
})

test_that("kmcompareplot model curve is the mean of the individual predicted curves", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype + karno, data = veteran)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  k <- kmcompareplot(fit, by = "celltype")
  expect_setequal(unique(k$group), levels(veteran$celltype))
  m <- k[k$group == "adeno" & k$curve == "model", ]
  idx <- veteran$celltype == "adeno"
  pp <- predict(fit, newdata = veteran[idx, ], times = m$time[-1])
  expect_equal(m$surv, c(1, unname(tapply(pp$est, pp$time, mean))), tolerance = 1e-12)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = veteran[idx, ])
  expect_equal(k$surv[k$group == "adeno" & k$curve == "km"], c(1, km$surv))

  expect_setequal(unique(kmcompareplot(fit, by = cut(veteran$karno, 3))$group),
                  levels(cut(veteran$karno, 3)))
  expect_error(kmcompareplot(fit, by = "karno"), "distinct values")
  expect_error(kmcompareplot(fit, by = "foo"), "must name a variable")
})

test_that("kmcompareplot model curve is close to KM for a well-specified model", {
  set.seed(11)
  n <- 4000
  x <- rnorm(n)
  z <- factor(sample(c("a", "b"), n, TRUE))
  t_ev <- (-log(runif(n)) / (0.1 * exp(1.2 * x + 0.7 * (z == "b"))))^(1 / 1.5)
  cens <- runif(n, 0, 15)
  d <- data.frame(time = pmin(t_ev, cens), status = as.numeric(t_ev <= cens), x = x, z = z)
  fit <- rpsurv(survival::Surv(time, status) ~ x + z, data = d)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  k <- kmcompareplot(fit, by = "z")
  for (g in c("a", "b")) {
    m <- k[k$group == g & k$curve == "model", ]
    km <- k[k$group == g & k$curve == "km", ]
    km_at <- stats::stepfun(km$time[-1], km$surv)(m$time)
    keep <- m$time < 10
    expect_lt(max(abs(m$surv[keep] - km_at[keep])), 0.04)
    # the curve at the mean covariate profile is a different estimand
    mean_prof <- predict(fit, newdata = data.frame(x = mean(d$x[d$z == g]), z = g),
                         times = m$time[keep][-1])$est
    expect_gt(max(abs(mean_prof - km_at[keep][-1])), 0.05)
  }
})

test_that("kmcompareplot handles left truncation through the risk set", {
  set.seed(7)
  n <- 4000
  x <- rnorm(n)
  t_ev <- (-log(runif(n)) / (0.1 * exp(0.8 * x)))^(1 / 1.5)
  e <- runif(n, 0, 3)
  keep <- t_ev > e
  cens <- e + runif(n, 0, 12)
  d <- data.frame(entry = e, time = pmin(t_ev, cens), status = as.numeric(t_ev <= cens),
                  x = x)[keep, ]
  fit <- rpsurv(survival::Surv(entry, time, status) ~ x, data = d)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  k <- kmcompareplot(fit)
  m <- k[k$curve == "model", ]
  km <- k[k$curve == "km", ]
  km_at <- stats::stepfun(km$time[-1], km$surv)(m$time)
  expect_lt(max(abs(m$surv - km_at)[m$time < 10]), 0.03)

  # split follow-up with constant covariates gives the same model curve
  bc <- brcancer
  half <- bc$rectime / 2
  long <- rbind(data.frame(start = 0, stop = half, status = 0, hormon = bc$hormon),
                data.frame(start = half, stop = bc$rectime, status = bc$censrec, hormon = bc$hormon))
  f_cp <- rpsurv(survival::Surv(start, stop, status) ~ hormon, data = long)
  f_rc <- rpsurv(survival::Surv(rectime, censrec) ~ hormon, data = bc)
  idx <- which(bc$hormon == 1)
  grid <- seq(0, max(bc$rectime), length.out = 201)
  expect_equal(rpsurv:::rp_average_survival(f_cp, which(long$hormon == 1), grid, TRUE),
               rpsurv:::rp_average_survival(f_rc, idx, grid, TRUE), tolerance = 1e-3)
})

test_that("coxsnellplot uses Nelson-Aalen and a square frame with the y = x line", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  cs <- coxsnellplot(fit)
  na <- survival::survfit(survival::Surv(residuals(fit), veteran$status) ~ 1, ctype = 1)
  expect_equal(cs$cumhaz, na$cumhaz)
  usr <- graphics::par("usr")
  expect_equal(usr[1:2], usr[3:4])
  expect_true(all(is.finite(cs$cumhaz)))
})

test_that("plot.tvecurve labels match what tvecurve() returns on each scale", {
  for (sc in c("hazard", "odds", "normal")) {
    fit <- rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran, tve = "karno",
                  scale = sc)
    tv <- tvecurve(fit, "karno", ref = 50, times = c(30, 100))
    expect_equal(rpsurv:::rp_tve_label(tv, TRUE), "Hazard ratio, karno = 51 vs 50")
    hr <- predict(fit, newdata = data.frame(trt = 0, karno = 51),
                  newdata0 = data.frame(trt = 0, karno = 50), times = c(30, 100), type = "hr")
    expect_equal(exp(tv$est), hr$est, tolerance = 1e-10)
    tl <- tvecurve(fit, "karno", type = "link")
    expect_match(rpsurv:::rp_tve_label(tl, FALSE),
                 switch(sc, hazard = "log cumulative hazard ratio", odds = "log cumulative odds ratio",
                        normal = "probit-scale effect"))
  }
})

test_that("plot(type = 'coxsnell' / 'km') dispatches to the diagnostic plots", {
  fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype, data = veteran, df = 3)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(plot(fit, type = "coxsnell"), coxsnellplot(fit))
  expect_identical(plot(fit, type = "km", by = "celltype"), kmcompareplot(fit, by = "celltype"))
  expect_identical(plot(fit, type = "km"), kmcompareplot(fit))
  expect_error(plot(fit, type = "coxsnell", newdata = veteran[1, ]), "does not apply")
  expect_error(plot(fit, type = "km", times = 1:10, ci = FALSE), "do not apply")
  expect_error(plot(fit, type = "coxsnell", by = "celltype"), "type = \"km\" only")
  expect_error(plot(fit, by = "celltype"), "type = \"km\" only")
})
