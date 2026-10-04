#' Non-linear effect component of a covariate
#'
#' Returns the fitted non-linear effect g(x) - g(ref) of a covariate named in
#' `nle`, on the linear-predictor scale, with a delta-method 95% confidence
#' interval. Because the spline enters the linear predictor additively and
#' does not depend on time, this is an exact linear combination of the
#' coefficients. On the hazard scale, and when `var` is not also in `tve`, it
#' is the log hazard ratio of `x` versus `ref` at every time. When `var` is
#' also in `tve`, it is the time-constant component of the effect, separate
#' from the time-varying component returned by [tvecurve()].
#'
#' @param object a fitted `"rpsurv"` object.
#' @param var name of a covariate in `object$nle`.
#' @param at numeric values of the covariate at which to evaluate g. Defaults
#'   to 50 points over the observed range.
#' @param ref reference value; `g(ref)` is subtracted. Defaults to the median
#'   of the covariate in the fitted data.
#' @param se.fit logical; add `se`, `lower`, and `upper`.
#' @return a data frame with columns `x` and `est`, and optionally `se`,
#'   `lower`, and `upper`.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' x <- runif(n, 0, 3)
#' dat <- data.frame(time = rexp(n, 0.2 * exp((x - 1.5)^2 / 2)),
#'                   status = 1, x = x)
#' fit <- rpsurv(survival::Surv(time, status) ~ x, data = dat, nle = "x")
#' head(nlecurve(fit, "x", ref = 1.5))
#' }
#' @export
nlecurve <- function(object, var, at = NULL, ref = NULL, se.fit = TRUE) {
  stopifnot(inherits(object, "rpsurv"), length(var) == 1L)
  if (is.null(object$nle) || !var %in% object$nle) {
    stop("`var` must be a covariate fitted with `nle`", call. = FALSE)
  }
  x_fit <- as.numeric(object$cov_data[[var]])
  if (is.null(at)) at <- seq(min(x_fit), max(x_fit), length.out = 50)
  if (is.null(ref)) ref <- stats::median(x_fit)
  k <- object$nle_knots[[var]]
  L <- rp_nle_basis(at, k) - rp_nle_basis(rep(ref, length(at)), k)
  cols <- paste0(var, ":s(x)", seq_len(ncol(L)))
  est <- as.numeric(L %*% object$coefficients[cols])
  out <- data.frame(x = at, est = est)
  if (se.fit) {
    se <- rp_delta_se(L, object$vcov[cols, cols, drop = FALSE])
    z <- stats::qnorm(0.975)
    out$se <- se
    out$lower <- est - z * se
    out$upper <- est + z * se
  }
  out
}

#' Time-varying effect component of a covariate
#'
#' Returns the fitted time-varying effect of a covariate named in `tve`, with
#' a delta-method 95% confidence interval.
#'
#' * `type = "hr"` (default): the log of the instantaneous hazard ratio at
#'   each time for `var = ref + 1` versus `var = ref`, as in
#'   [predict.rpsurv()] with `type = "hr"`. Other covariates are set to 0;
#'   on the hazard scale they cancel unless they also have a time-varying
#'   effect. If `var` is also in `nle`, its non-linear component is held
#'   fixed, so the result is the time-varying component only (see
#'   [nlecurve()] for the other one).
#' * `type = "link"`: the coefficient function beta(t) on the linear-predictor
#'   (log cumulative hazard, log cumulative odds, or probit) scale, the main
#'   effect plus its spline in log time; when `var` is also in `nle` the main
#'   effect belongs to g, so beta(t) is the spline part only.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param var name of a covariate in `object$tve`.
#' @param times times at which to evaluate the effect. Defaults to 50 points
#'   over the observed event times.
#' @param type `"hr"` or `"link"`.
#' @param ref reference value of `var` for `type = "hr"`. Default 0.
#' @param se.fit logical; add `se`, `lower`, and `upper`.
#' @return a data frame with columns `time` and `est`, and optionally `se`,
#'   `lower`, and `upper`.
#' @examples
#' \donttest{
#' fit <- rpsurv(survival::Surv(time, status) ~ trt + karno, data = veteran,
#'               tve = "karno", nle = "karno")
#' tvecurve(fit, "karno", times = c(30, 90, 180), ref = 50)
#' }
#' @export
tvecurve <- function(object, var, times = NULL, type = c("hr", "link"), ref = 0,
                     se.fit = TRUE) {
  stopifnot(inherits(object, "rpsurv"), length(var) == 1L)
  type <- match.arg(type)
  if (is.null(object$tve) || !var %in% object$tve) {
    stop("`var` must be a covariate fitted with `tve`", call. = FALSE)
  }
  if (is.null(times)) {
    ev <- object$time[object$status == 1]
    times <- seq(min(ev), max(ev), length.out = 50)
  }
  log_times <- log(times)
  b <- object$coefficients
  V <- object$vcov
  in_nle <- !is.null(object$nle) && var %in% object$nle
  nle_cols <- if (in_nle) grep(paste0("^", var, ":s\\(x\\)"), names(b)) else integer(0)

  if (type == "link") {
    tcols <- grep(paste0("^", var, ":s\\(logt\\)"), names(b))
    L <- matrix(0, length(times), length(b))
    L[, tcols] <- rcs_basis(log_times, object$tve_knots[[var]])
    if (!in_nle) L[, match(var, names(b))] <- 1
    est <- as.numeric(L %*% b)
    se <- rp_delta_se(L, V)
  } else {
    profile <- function(value) {
      cov <- as.data.frame(lapply(object$cov_data, function(v) 0))
      cov[[var]] <- value
      d <- rp_predict_design(object, cov, log_times)
      if (in_nle) {
        d$X[, nle_cols] <- 0
        d$eta <- as.numeric(d$X %*% b)
        d$S <- rp_inv_link(object$scale)(d$eta)
      }
      d
    }
    d1 <- profile(ref + 1)
    d0 <- profile(ref)
    est <- log(rp_hazard_from_design(d1, times, object$scale) /
                 rp_hazard_from_design(d0, times, object$scale))
    grad <- rp_grad_loghaz(d1, object$scale) - rp_grad_loghaz(d0, object$scale)
    if (in_nle) grad[, nle_cols] <- 0
    se <- rp_delta_se(grad, V)
  }
  out <- data.frame(time = times, est = est)
  if (se.fit) {
    z <- stats::qnorm(0.975)
    out$se <- se
    out$lower <- est - z * se
    out$upper <- est + z * se
  }
  out
}
