rp_effect_label <- function(scale) {
  switch(scale, hazard = "Hazard ratio", odds = "Odds ratio", normal = "Probit effect")
}

rp_plot_effect <- function(xv, d, exp, ci, log, col, ref.line, add, xlab, ylab, ...) {
  tr <- if (exp) base::exp else identity
  est <- tr(d$est)
  has_ci <- ci && !is.null(d$lower)
  lo <- if (has_ci) tr(d$lower) else est
  hi <- if (has_ci) tr(d$upper) else est
  if (!add) {
    yr <- range(c(lo, hi, if (ref.line) as.numeric(exp)), finite = TRUE)
    graphics::plot(range(xv), yr, type = "n", log = log, xlab = xlab, ylab = ylab, ...)
    if (ref.line) graphics::abline(h = as.numeric(exp), lty = 2, col = "grey50")
  }
  if (has_ci) {
    ok <- is.finite(lo) & is.finite(hi)
    graphics::polygon(c(xv[ok], rev(xv[ok])), c(lo[ok], rev(hi[ok])),
                      col = grDevices::adjustcolor(col, alpha.f = 0.15), border = NA)
  }
  graphics::lines(xv, est, col = col, lwd = 2)
  invisible(d)
}

#' Plot a time-varying effect
#'
#' Draws the output of [tvecurve()]: the effect against time, with its 95%
#' confidence band and a reference line at no effect.
#'
#' @param x a `"tvecurve"` object returned by [tvecurve()].
#' @param exp logical; exponentiate the estimate and band (hazard ratio
#'   instead of log hazard ratio). Defaults to `TRUE` for `type = "hr"` and
#'   `FALSE` for `type = "link"`.
#' @param ci logical; draw the confidence band when available.
#' @param log passed to [graphics::plot()]; `"y"` (log y axis) when `exp`
#'   is `TRUE`, `""` otherwise.
#' @param col colour of the curve and band.
#' @param ref.line logical; draw a dashed line at 1 (`exp = TRUE`) or 0.
#' @param add logical; add to the current plot instead of starting a new one.
#' @param xlab,ylab axis labels.
#' @param ... further arguments passed to [graphics::plot()].
#' @return Called for its side effect. Invisibly returns `x`.
#' @examples
#' \donttest{
#' fit <- rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran,
#'               tve = "karno")
#' plot(tvecurve(fit, "karno", ref = 50))
#' }
#' @export
plot.tvecurve <- function(x, exp = attr(x, "type") != "link", ci = TRUE,
                          log = if (exp) "y" else "", col = 1, ref.line = TRUE,
                          add = FALSE, xlab = "Time", ylab = NULL, ...) {
  if (is.null(ylab)) {
    var <- attr(x, "var")
    ylab <- if (identical(attr(x, "type"), "link")) {
      paste0("beta(t) for ", var)
    } else {
      paste0(if (exp) "Hazard ratio" else "Log hazard ratio", " per unit of ", var)
    }
  }
  rp_plot_effect(x$time, x, exp, ci, log, col, ref.line, add, xlab, ylab, ...)
  invisible(x)
}

#' Plot a non-linear effect
#'
#' Draws the output of [nlecurve()]: g(x) - g(ref) against the covariate,
#' with its 95% confidence band and a reference line at no effect.
#'
#' @inheritParams plot.tvecurve
#' @param x a `"nlecurve"` object returned by [nlecurve()].
#' @param exp logical; exponentiate the estimate and band. Defaults to
#'   `TRUE` on the hazard and odds scales, `FALSE` on the normal scale.
#' @param rug logical; mark the reference value on the x axis.
#' @return Called for its side effect. Invisibly returns `x`.
#' @examples
#' \donttest{
#' fit <- rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran,
#'               nle = "karno")
#' plot(nlecurve(fit, "karno", ref = 50))
#' }
#' @export
plot.nlecurve <- function(x, exp = !identical(attr(x, "scale"), "normal"), ci = TRUE,
                          log = if (exp) "y" else "", col = 1, ref.line = TRUE,
                          add = FALSE, xlab = NULL, ylab = NULL, rug = TRUE, ...) {
  var <- attr(x, "var")
  ref <- attr(x, "ref")
  if (is.null(xlab)) xlab <- var
  if (is.null(ylab)) {
    lab <- rp_effect_label(attr(x, "scale"))
    if (!exp && lab != "Probit effect") lab <- paste("Log", tolower(lab))
    ylab <- paste0(lab, " vs ", var, " = ", signif(ref, 3))
  }
  rp_plot_effect(x$x, x, exp, ci, log, col, ref.line, add, xlab, ylab, ...)
  if (rug && !add) graphics::rug(ref, col = col, lwd = 2)
  invisible(x)
}
