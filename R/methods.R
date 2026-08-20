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

#' @importFrom stats AIC
#' @export
AIC.rpsurv <- function(object, ..., k = 2) -2 * object$loglik + k * object$df

#' @importFrom stats BIC
#' @export
BIC.rpsurv <- function(object, ...) -2 * object$loglik + log(object$nevent) * object$df

#' @importFrom stats confint
#' @export
confint.rpsurv <- function(object, parm, level = 0.95, ...) {
  se <- sqrt(diag(object$vcov))
  b <- object$coefficients
  if (missing(parm)) parm <- names(b)
  z <- stats::qnorm(1 - (1 - level) / 2)
  lo <- b[parm] - z * se[parm]
  hi <- b[parm] + z * se[parm]
  out <- cbind(lo, hi)
  colnames(out) <- sprintf(c("%.1f %%", "%.1f %%"), c((1 - level) / 2 * 100, (1 + level) / 2 * 100))
  out
}

is_baseline_term <- function(object) {
  seq_along(object$coefficients) <= object$p_base
}

#' @export
print.rpsurv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Royston-Parmar flexible parametric survival model\n")
  cat("Call:\n")
  print(x$call)
  cat("\nScale:", x$scale, "  Baseline df:", x$p_base - 1L,
      if (!is.null(x$tve)) paste0("  Time-varying effect: ", paste(x$tve, collapse = ", ")) else "", "\n")
  if (isTRUE(x$counting)) cat("Data: counting-process (left truncation / time-varying covariates)\n")
  cat("n =", x$n, ", number of events =", x$nevent, "\n")
  cat("Log-likelihood =", format(x$loglik, digits = digits),
      "  AIC =", format(AIC(x), digits = digits), "\n")
  invisible(x)
}

#' @export
summary.rpsurv <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
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
      counting = object$counting
    ),
    class = "summary.rpsurv"
  )
}

#' @export
print.summary.rpsurv <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  cat("Royston-Parmar flexible parametric survival model\n")
  cat("Call:\n")
  print(x$call)

  hr_label <- switch(x$scale, hazard = "exp(coef) [HR]", odds = "exp(coef) [OR]", "exp(coef)")

  cov_tab <- x$coefficients[!x$baseline, , drop = FALSE]
  if (nrow(cov_tab)) {
    colnames(cov_tab)[2] <- hr_label
    cat("\nCovariate effects (", x$scale, " scale):\n", sep = "")
    stats::printCoefmat(cov_tab, digits = digits, signif.stars = signif.stars,
                         P.values = TRUE, has.Pvalue = TRUE)
  } else {
    cat("\n(no covariates in the model)\n")
  }

  base_tab <- x$coefficients[x$baseline, , drop = FALSE]
  cat("\nBaseline spline terms:\n")
  stats::printCoefmat(base_tab, digits = digits, signif.stars = FALSE,
                       P.values = TRUE, has.Pvalue = TRUE)

  cat("\nn =", x$n, ", number of events =", x$nevent, ", parameters =", x$df, "\n")
  cat("Log-likelihood =", format(x$loglik, digits = digits),
      "  AIC =", format(x$aic, digits = digits),
      "  BIC =", format(x$bic, digits = digits), "\n")
  invisible(x)
}
