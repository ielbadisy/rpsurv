#' @export
coef.rpsurv <- function(object, ...) object$coefficients

#' @export
vcov.rpsurv <- function(object, ...) object$vcov

#' @export
logLik.rpsurv <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$df
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

# One value per model, or a data frame with df and the criterion (one row per model, named as in the
# call) when several are given, as stats::AIC() and stats::BIC(). Other model classes use their own method.
rp_ic <- function(ms, nms, pen, ic, f) {
  v <- vapply(ms, function(m) if (inherits(m, "rpsurv")) -2 * m$loglik + pen(m) * m$df else f(m), 0)
  if (length(ms) == 1L) return(v)
  df <- vapply(ms, function(m) if (inherits(m, "rpsurv")) m$df else attr(stats::logLik(m), "df"), 0)
  stats::setNames(data.frame(df, v, row.names = nms), c("df", ic))
}

#' @importFrom stats AIC
#' @export
AIC.rpsurv <- function(object, ..., k = 2) {
  nms <- vapply(as.list(match.call(expand.dots = TRUE)[-1L]), deparse1, "")
  rp_ic(list(object, ...), nms[names(nms) != "k"], function(m) k, "AIC", function(m) stats::AIC(m, k = k))
}

# The BIC penalty uses the number of events, the effective sample size of censored data (Volinsky and
# Raftery 2000).
#' @importFrom stats BIC
#' @export
BIC.rpsurv <- function(object, ...) {
  nms <- vapply(as.list(match.call(expand.dots = TRUE)[-1L]), deparse1, "")
  rp_ic(list(object, ...), nms, function(m) log(m$nevent), "BIC", stats::BIC)
}

# Standard errors from the vcov diagonal. A negative or non-finite variance
# (non-positive-definite vcov, already reported by rpsurv()) gives NA plus
# one clear warning instead of R's generic "NaNs produced".
rp_coef_se <- function(object) {
  v <- diag(object$vcov)
  bad <- !is.finite(v) | v < 0
  if (any(bad)) {
    warning("Standard errors are NA for ", paste(names(object$coefficients)[bad], collapse = ", "),
            ": the variance-covariance matrix is not positive definite.", call. = FALSE)
    v[bad] <- NA_real_
  }
  sqrt(v)
}

#' @importFrom stats confint
#' @export
confint.rpsurv <- function(object, parm, level = 0.95, ...) {
  se <- rp_coef_se(object)
  b <- object$coefficients
  if (missing(parm)) parm <- names(b)
  z <- stats::qnorm(1 - (1 - level) / 2)
  lo <- b[parm] - z * se[parm]
  hi <- b[parm] + z * se[parm]
  out <- cbind(lo, hi)
  colnames(out) <- sprintf(c("%.1f %%", "%.1f %%"), c((1 - level) / 2 * 100, (1 + level) / 2 * 100))
  out
}

# "(k observations deleted due to missingness)", as in print.coxph()
rp_na_note <- function(na_action) {
  k <- length(na_action)
  if (k) sprintf(" (%d observation%s deleted due to missingness)", k, if (k > 1L) "s" else "") else ""
}

# Log-likelihood, AIC and BIC with two decimals, so models can be compared
# and LR tests done by hand; format() keeps the session's digits and scipen.
rp_fmt_fit <- function(v) format(round(v, 2), nsmall = 2)

rp_title <- function(x) {
  if (is.null(x$cause)) return("Royston-Parmar flexible parametric survival model")
  paste0("Royston-Parmar flexible parametric Fine-Gray model (cause ", x$cause, ", robust SE)")
}

is_baseline_term <- function(object) {
  seq_along(object$coefficients) <= object$p_base
}

#' @export
print.rpsurv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(rp_title(x), "\n", sep = "")
  cat("Call:\n")
  print(x$call)
  cat("\nScale: ", x$scale, "   Baseline df: ", x$p_base - 1L,
      if (!is.null(x$tve)) paste0("   Time-varying effect: ", paste(x$tve, collapse = ", ")),
      if (!is.null(x$nle)) paste0("   Non-linear effect: ", paste(x$nle, collapse = ", ")),
      "\n", sep = "")
  if (isTRUE(x$counting)) cat("Data: counting-process (left truncation / time-varying covariates)\n")
  cat("n = ", x$n, ", number of events = ", x$nevent, rp_na_note(x$na.action), "\n", sep = "")
  cat("Log-likelihood = ", rp_fmt_fit(x$loglik), ", AIC = ", rp_fmt_fit(AIC(x)), "\n", sep = "")
  invisible(x)
}

#' @export
summary.rpsurv <- function(object, ...) {
  se <- rp_coef_se(object)
  b <- object$coefficients
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(coef = b, `exp(coef)` = exp(b), `se(coef)` = se, z = z, `Pr(>|z|)` = p)
  baseline <- is_baseline_term(object)

  structure(
    list(
      call = object$call,
      coefficients = tab,
      baseline = baseline,
      scale = object$scale,
      loglik = object$loglik,
      n = object$n,
      nevent = object$nevent,
      df = object$df,
      aic = AIC(object),
      bic = BIC(object),
      tve = object$tve,
      nle = object$nle,
      counting = object$counting,
      na.action = object$na.action,
      cause = object$cause
    ),
    class = "summary.rpsurv"
  )
}

#' @export
print.summary.rpsurv <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  cat(rp_title(x), "\n", sep = "")
  cat("Call:\n")
  print(x$call)

  hr_label <- switch(x$scale, hazard = if (is.null(x$cause)) "exp(coef) [HR]" else "exp(coef) [sHR]", odds = "exp(coef) [OR]", "exp(coef)")

  cov_tab <- x$coefficients[!x$baseline, , drop = FALSE]
  if (nrow(cov_tab)) {
    colnames(cov_tab)[2] <- hr_label
    cat("\nCovariate effects (", x$scale, " scale):\n", sep = "")
    stats::printCoefmat(cov_tab, digits = digits, signif.stars = signif.stars,
                         P.values = TRUE, has.Pvalue = TRUE)
  } else {
    cat("\n(no covariates in the model)\n")
  }

  # exp() of a spline coefficient has no interpretation, and a tiny value
  # such as exp(intercept) would push the whole table into scientific notation
  base_tab <- x$coefficients[x$baseline, colnames(x$coefficients) != "exp(coef)", drop = FALSE]
  cat("\nBaseline spline terms:\n")
  stats::printCoefmat(base_tab, digits = digits, signif.stars = FALSE,
                       P.values = TRUE, has.Pvalue = TRUE)

  cat("\nn = ", x$n, ", number of events = ", x$nevent, ", parameters = ", x$df,
      rp_na_note(x$na.action), "\n", sep = "")
  cat("Log-likelihood = ", rp_fmt_fit(x$loglik), ", AIC = ", rp_fmt_fit(x$aic),
      ", BIC = ", rp_fmt_fit(x$bic), "\n", sep = "")
  invisible(x)
}
