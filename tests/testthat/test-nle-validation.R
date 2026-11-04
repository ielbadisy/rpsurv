# Validation of non-linear effects (`nle`) against rstpm2::stpm2. A natural
# spline ns(x, knots, Boundary.knots) with rpsurv's nle knots spans the same
# model space as rpsurv's restricted cubic spline in x, so the two fits must
# agree exactly once rpsurv's coefficients are mapped into stpm2's basis.

ns_term <- function(fit, v) {
  k <- fit$nle_knots[[v]]
  sprintf("splines::ns(%s, knots = c(%s), Boundary.knots = c(%s))", v,
          paste(format(k[-c(1, length(k))], digits = 17), collapse = ", "),
          paste(format(range(k), digits = 17), collapse = ", "))
}

rp_in_stpm2_basis <- function(fit, ref) {
  M <- qr.solve(ref@args$X, fit$X)
  list(coef = drop(M %*% coef(fit)), vcov = M %*% vcov(fit) %*% t(M))
}

stpm2_twin <- function(fit, lhs, terms, data, link) {
  rhs <- vapply(terms, function(v) if (v %in% fit$nle) ns_term(fit, v) else v, "")
  f <- stats::as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
  tvc <- if (length(fit$tve)) {
    stats::setNames(lapply(fit$tve, function(v) length(fit$tve_knots[[v]]) - 1L), fit$tve)
  }
  if (is.null(tvc)) {
    rstpm2::stpm2(f, data = data, df = 4, link.type = link)
  } else {
    rstpm2::stpm2(f, data = data, df = 4, tvc = tvc, link.type = link)
  }
}

test_that("nle fits match rstpm2::stpm2 in log-likelihood, coefficients and full vcov", {
  skip_if_not_installed("rstpm2")
  skip_on_cran()
  suppressPackageStartupMessages(library(rstpm2)) # stpm2() calls internal gsm() unqualified
  d <- get(data("brcancer", package = "rstpm2", envir = environment()))
  d$hormon <- as.numeric(d$hormon)
  links <- c(hazard = "PH", odds = "PO", normal = "probit")
  for (sc in names(links)) {
    cases <- list(
      list(rpsurv(survival::Surv(time, status) ~ trt + karno + age, data = veteran,
                  nle = c("age", "karno"), scale = sc),
           "survival::Surv(time, status)", c("trt", "karno", "age"), veteran),
      list(rpsurv(survival::Surv(time, status) ~ trt + karno + age, data = veteran,
                  nle = c("age", "karno"), tve = "karno", scale = sc),
           "survival::Surv(time, status)", c("trt", "karno", "age"), veteran),
      list(rpsurv(survival::Surv(rectime, censrec) ~ hormon + x1 + x6, data = d,
                  nle = c("x1", "x6"), scale = sc),
           "survival::Surv(rectime, censrec)", c("hormon", "x1", "x6"), d),
      list(rpsurv(survival::Surv(rectime, censrec) ~ hormon + x1 + x6, data = d,
                  nle = c("x1", "x6"), tve = "x1", tve.df = 2, scale = sc),
           "survival::Surv(rectime, censrec)", c("hormon", "x1", "x6"), d)
    )
    for (cs in cases) {
      fit <- cs[[1]]
      ref <- stpm2_twin(fit, cs[[2]], cs[[3]], cs[[4]], links[[sc]])
      M <- qr.solve(ref@args$X, fit$X)
      expect_lt(max(abs(ref@args$X %*% M - fit$X)), 1e-8) # same model space
      mapped <- rp_in_stpm2_basis(fit, ref)
      expect_equal(fit$loglik, -ref@min, tolerance = 1e-7)
      expect_equal(unname(mapped$coef), unname(coef(ref)), tolerance = 1e-3)
      expect_equal(unname(mapped$vcov), unname(vcov(ref)), tolerance = 1e-3)
    }
  }
})

test_that("nlecurve and tvecurve agree with a numerical delta method", {
  skip_if_not_installed("numDeriv")
  withcoef <- function(fit, b) {
    fit$coefficients[] <- b
    fit
  }
  num_se <- function(fit, fn) {
    G <- numDeriv::jacobian(fn, coef(fit))
    sqrt(rowSums((G %*% vcov(fit)) * G))
  }
  at <- c(40, 50, 70, 80)
  tt <- c(10, 50, 100, 300)
  for (sc in c("hazard", "odds", "normal")) {
    fit <- rpsurv(survival::Surv(time, status) ~ trt + karno + age, data = veteran,
                  nle = "age", tve = "karno", scale = sc)
    g <- nlecurve(fit, "age", at = at, ref = 62)
    link_diff <- function(b) {
      f <- withcoef(fit, b)
      predict(f, newdata = data.frame(trt = 1, karno = 60, age = at), times = 50, type = "link")$est -
        predict(f, newdata = data.frame(trt = 1, karno = 60, age = 62), times = 50, type = "link")$est
    }
    expect_equal(g$est, link_diff(coef(fit)), tolerance = 1e-10)
    expect_equal(g$se, num_se(fit, link_diff), tolerance = 1e-6)

    tv <- tvecurve(fit, "karno", times = tt, ref = 50)
    log_hr <- function(b) {
      f <- withcoef(fit, b)
      log(predict(f, newdata = data.frame(trt = 0, karno = 51, age = 0),
                  newdata0 = data.frame(trt = 0, karno = 50, age = 0), times = tt, type = "hr")$est)
    }
    expect_equal(tv$se, num_se(fit, log_hr), tolerance = 1e-6)
  }
})
