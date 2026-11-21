mgus <- function() {
  d <- survival::mgus2
  d$etime <- with(d, ifelse(pstat == 1, ptime, futime))
  d$event <- factor(with(d, ifelse(pstat == 1, 1, 2 * death)), 0:2, c("censor", "pcm", "death"))
  d
}

test_that("integer case weights give the same fit as replicated rows", {
  set.seed(1)
  d <- veteran
  w <- sample(1:3, nrow(d), replace = TRUE)
  fw <- rpsurv(survival::Surv(time, status) ~ karno + trt, data = d, weights = w)
  # same knots: the default ones are centiles of the (unweighted) event times
  fr <- rpsurv(survival::Surv(time, status) ~ karno + trt, data = d[rep(seq_len(nrow(d)), w), ],
               knots = fw$knots)
  expect_equal(coef(fw), coef(fr), tolerance = 1e-5)
  expect_equal(fw$loglik, fr$loglik, tolerance = 1e-6)
  expect_equal(vcov(fw), vcov(fr), tolerance = 1e-3)
  d$w <- w
  expect_equal(coef(rpsurv(survival::Surv(time, status) ~ karno + trt, data = d, weights = w)), coef(fw))
})

test_that("row scores add up to the C++ gradient on every scale", {
  d <- veteran
  d$start <- pmin(d$time / 2, 5)
  for (sc in c("hazard", "odds", "normal")) {
    f <- rpsurv(survival::Surv(start, time, status) ~ karno, data = d, scale = sc, df = 2)
    f$coefficients <- f$coefficients + 0.01
    w <- runif(nrow(d))
    g <- rp_negloglik_grad_cpp(f$coefficients, f$X, f$dX, f$Xentry, f$has_entry, f$log_time,
                               as.numeric(f$status), rp_scale_code(sc), w)$gradient
    expect_equal(colSums(w * rp_scores(f)), -g, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("a Cox fit on the IPCW expansion gives the Fine-Gray estimates of crr and finegray", {
  d <- mgus()
  e <- rp_fgexpand(d$etime, as.integer(d$event) - 1L, 1L)
  fg <- cbind(d[e$id, ], e)
  b <- coef(survival::coxph(survival::Surv(start, stop, status) ~ sex + age, weights = w, data = fg))
  ref <- survival::finegray(survival::Surv(etime, event) ~ sex + age, data = d, etype = "pcm")
  bf <- coef(survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ sex + age, weights = fgwt,
                             data = ref))
  expect_equal(b, bf, tolerance = 5e-3)
  skip_if_not_installed("cmprsk")
  cc <- cmprsk::crr(d$etime, as.integer(d$event) - 1L, cbind(d$sex == "M", d$age), failcode = 1, cencode = 0)
  expect_equal(unname(b), unname(cc$coef), tolerance = 1e-3)
})

test_that("subdistribution HRs and robust SEs agree with cmprsk::crr", {
  skip_if_not_installed("cmprsk")
  d <- mgus()
  fit <- rpcr(survival::Surv(etime, event) ~ sex + age, data = d, cause = "pcm")
  cc <- cmprsk::crr(d$etime, as.integer(d$event) - 1L, cbind(sexM = d$sex == "M", age = d$age),
                    failcode = 1, cencode = 0)
  b <- coef(fit)[c("sexM", "age")]
  expect_equal(unname(exp(b)), unname(exp(cc$coef)), tolerance = 0.01)
  expect_equal(unname(sqrt(diag(vcov(fit)))[c("sexM", "age")]), sqrt(diag(cc$var)), tolerance = 0.02)
  expect_s3_class(fit, c("rpcr", "rpsurv"))
  expect_equal(fit$n, nrow(d))
  expect_equal(fit$nevent, sum(d$event == "pcm"))
})

test_that("the cumulative incidence agrees with Aalen-Johansen", {
  d <- mgus()
  tt <- c(12, 60, 120, 240)
  aj <- summary(survival::survfit(survival::Surv(etime, event) ~ 1, data = d), times = tt)
  for (k in c("pcm", "death")) {
    p <- predict(rpcr(survival::Surv(etime, event) ~ 1, data = d, cause = k), data.frame(x = 1),
                 times = tt, se.fit = TRUE)
    expect_lt(max(abs(p$est - aj$pstate[, match(k, aj$states)])), 0.01)
    expect_true(all(p$lower < p$est & p$est < p$upper & p$lower >= 0 & p$upper <= 1))
  }
})

test_that("predict types are consistent", {
  d <- mgus()
  fit <- rpcr(survival::Surv(etime, event) ~ sex + age, data = d, cause = "pcm")
  n1 <- data.frame(sex = "M", age = 70); n0 <- data.frame(sex = "F", age = 70)
  tt <- c(60, 240)
  c1 <- predict(fit, n1, times = tt); c0 <- predict(fit, n0, times = tt)
  dd <- predict(fit, n1, n0, times = tt, type = "cifdiff", se.fit = TRUE)
  expect_equal(dd$est, c1$est - c0$est)
  expect_true(all(dd$lower < dd$est & dd$est < dd$upper))
  s <- predict(fit, n1, n0, times = tt, type = "shr", se.fit = TRUE)
  expect_equal(s$est, rep(unname(exp(coef(fit)["sexM"])), 2))
  expect_equal(s$upper, rep(unname(exp(confint(fit)["sexM", 2])), 2))
  expect_equal(nrow(predict(fit, times = 60)), nrow(d))
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_s3_class(plot(fit, newdata = rbind(n0, n1)), "data.frame")
  expect_s3_class(plot(fit, type = "shr", newdata = n1, newdata0 = n0), "data.frame")
  expect_output(summary(fit) |> print(), "sHR")
})

test_that("AIC and BIC of several models give one row per model, as stats::AIC", {
  f1 <- rpsurv(survival::Surv(time, status) ~ karno, data = veteran)
  f2 <- rpsurv(survival::Surv(time, status) ~ karno + trt, data = veteran)
  a <- AIC(f1, f2)
  expect_equal(dim(a), c(2L, 2L))
  expect_equal(rownames(a), c("f1", "f2"))
  expect_equal(a$df, c(f1$df, f2$df))
  expect_equal(a$AIC, c(AIC(f1), AIC(f2)))
  expect_equal(AIC(f1), -2 * f1$loglik + 2 * f1$df)
  expect_equal(AIC(f1, f2, k = 3)$AIC, c(AIC(f1, k = 3), AIC(f2, k = 3)))
  b <- BIC(f1, f2)
  expect_equal(names(b), c("df", "BIC"))
  expect_equal(b$BIC, c(BIC(f1), BIC(f2)))
  expect_equal(BIC(f1), -2 * f1$loglik + log(f1$nevent) * f1$df)
})

test_that("plot(type = 'cifdiff') draws the CIF difference; delta-method band matches simulation", {
  d <- mgus()
  fit <- rpcr(survival::Surv(etime, event) ~ sex, data = d, cause = "pcm")
  n1 <- data.frame(sex = "M"); n0 <- data.frame(sex = "F"); tt <- c(60, 120, 240)
  pdf(NULL); dev.control("enable"); on.exit(dev.off())
  p <- plot(fit, type = "cifdiff", newdata = n1, newdata0 = n0, times = tt)
  dd <- predict(fit, n1, n0, times = tt, type = "cifdiff", se.fit = TRUE)
  expect_equal(p, dd)
  expect_true(all(c("CIF difference (cause pcm)", "grey50") %in%
                    unlist(lapply(recordPlot()[[1]], function(e) Filter(is.character, e[[2]])))))
  expect_error(plot(fit, type = "cifdiff", newdata = n1), "newdata0")
  se <- (dd$upper - dd$lower) / (2 * qnorm(0.975))
  # delta method vs simulation from the asymptotic normal law of the coefficients
  set.seed(1)
  B <- sweep(matrix(rnorm(2000 * length(coef(fit))), 2000) %*% chol(vcov(fit)), 2, coef(fit), "+")
  sim <- apply(B, 1, function(b) {
    f <- fit; f$coefficients <- b
    predict(f, n1, n0, times = tt, type = "cifdiff")$est
  })
  expect_equal(se, apply(sim, 1, sd), tolerance = 0.06)
})

test_that("CIF difference agrees with the true difference and Aalen-Johansen under a Fine-Gray model", {
  # Fine (1999) generator: F1(t | x) = 1 - {1 - p (1 - exp(-t))}^exp(b x), known exactly
  F1 <- function(t, x, p = 0.4, b = 0.5) 1 - (1 - p * (1 - exp(-t)))^exp(b * x)
  set.seed(2); n <- 1500; p <- 0.4; b <- 0.5
  x <- rbinom(n, 1, 0.5); E <- exp(b * x); c1 <- runif(n) < 1 - (1 - p)^E
  t1 <- -log(1 - (1 - (1 - runif(n) * (1 - (1 - p)^E))^(1 / E)) / p)
  tt <- ifelse(c1, t1, rexp(n, exp(-0.5 * x))); cz <- runif(n, 0, 4)
  d <- data.frame(x, time = pmin(tt, cz), status = factor(ifelse(cz < tt, 0, ifelse(c1, 1, 2)), 0:2))
  fit <- rpcr(survival::Surv(time, status) ~ x, data = d, cause = "1")
  tg <- c(0.5, 1, 2)
  dd <- predict(fit, data.frame(x = 1), data.frame(x = 0), times = tg, type = "cifdiff", se.fit = TRUE)
  se <- (dd$upper - dd$lower) / (2 * qnorm(0.975))
  aj <- summary(survival::survfit(survival::Surv(time, status) ~ x, data = d), times = tg)
  a1 <- aj$strata == "x=1"; a0 <- aj$strata == "x=0"
  expect_lt(max(abs(dd$est - (aj$pstate[a1, 2] - aj$pstate[a0, 2]))), 0.01)
  expect_true(all(abs(dd$est - (F1(tg, 1) - F1(tg, 0))) < 2 * se))
  # the parametric model is at least as precise as the nonparametric estimate, of the same order
  ajse <- sqrt(aj$std.err[a1, 2]^2 + aj$std.err[a0, 2]^2)
  expect_true(all(se < 1.05 * ajse & se > 0.5 * ajse))
})
