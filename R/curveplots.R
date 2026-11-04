rp_effect_label <- function(scale) {
  switch(scale, hazard = "Hazard ratio", odds = "Odds ratio", normal = "Probit effect")
}

# Default y label of plot.tvecurve(): what tvecurve() returned, on which scale.
#' @keywords internal
rp_tve_label <- function(x, exp) {
  var <- attr(x, "var")
  scale <- attr(x, "scale")
  if (is.null(scale)) scale <- "hazard"
  if (identical(attr(x, "type"), "link")) {
    lab <- switch(scale,
      hazard = c("Log cumulative hazard ratio", "Cumulative hazard ratio"),
      odds = c("Log cumulative odds ratio", "Cumulative odds ratio"),
      normal = c("Probit-scale effect", "exp(probit-scale effect)")
    )[1L + exp]
    paste0("beta(t): ", tolower(substr(lab, 1, 1)), substring(lab, 2), " per unit of ", var)
  } else {
    ref <- attr(x, "ref")
    if (is.null(ref)) ref <- 0
    paste0(if (exp) "Hazard ratio" else "Log hazard ratio", ", ", var, " = ",
           signif(ref + 1, 4), " vs ", signif(ref, 4))
  }
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
#' @param exp logical; exponentiate the estimate and band. Defaults to
#'   `TRUE` for `type = "hr"` (a hazard ratio on every model scale, since
#'   [tvecurve()] computes it from the fitted hazards) and `FALSE` for
#'   `type = "link"`.
#' @details The default y label states what is drawn. For `type = "hr"` it
#'   is the hazard ratio of `var = ref + 1` versus `var = ref` (other
#'   covariates at 0). It is not a ratio "per unit" of `var`: with a
#'   time-varying effect the derivative of the linear predictor depends on
#'   `var`, so the hazard ratio changes with `ref`, and off the hazard scale
#'   it also depends on the other covariates. For `type = "link"` it is beta(t), the effect per unit of `var`
#'   on the model scale: a log cumulative hazard ratio (`scale = "hazard"`),
#'   a log cumulative odds ratio (`"odds"`), or a shift of the probit of the
#'   failure probability (`"normal"`); exponentiated, the first two are a
#'   cumulative hazard ratio and a cumulative odds ratio.
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
  if (is.null(ylab)) ylab <- rp_tve_label(x, exp)
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
