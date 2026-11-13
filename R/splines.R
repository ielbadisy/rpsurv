#' Restricted cubic spline basis (Royston-Parmar parameterisation)
#'
#' Builds the basis used by Royston & Parmar (2002) flexible parametric
#' models and by [rpsurv()] for the baseline log cumulative hazard, the
#' time-varying effects (`tve`) and the non-linear effects (`nle`): a
#' linear term plus one term per interior knot, each term being a natural
#' (restricted) cubic spline component that is linear beyond the boundary
#' knots. For knots \eqn{k_{min} < k_1 < \dots < k_m < k_{max}}, the column
#' for interior knot \eqn{k_j} is
#' \deqn{(x - k_j)_+^3 - \lambda_j (x - k_{min})_+^3 - (1 - \lambda_j)(x - k_{max})_+^3,
#'   \quad \lambda_j = (k_{max} - k_j) / (k_{max} - k_{min}).}
#'
#' Together with an intercept, the basis spans the same space as
#' `rstpm2::nsx(x, knots = <interior>, Boundary.knots = <boundary>)` and as
#' `rms::rcs(x, knots)` with the same knots; only the parameterisation
#' differs, so the coefficients differ but the fitted curves are identical.
#'
#' @param x numeric vector (typically log time) at which to evaluate the basis.
#' @param knots numeric vector of all knots: the two boundary knots and the
#'   interior knots, for example from [defaultknots()]. It is sorted
#'   internally.
#' @param derivative if `TRUE`, return the derivative of each basis column
#'   with respect to `x` instead of the basis itself.
#' @return a numeric matrix with `length(x)` rows and `length(knots) - 1`
#'   columns, the first being `x` itself (or 1 for the derivative).
#' @seealso [defaultknots()] for the default knot placement.
#' @examples
#' lt <- log(c(5, 10, 20, 40, 80, 160))
#' k <- defaultknots(lt, df = 3)
#' rcsbasis(lt, k)
#' rcsbasis(lt, k, derivative = TRUE)
#' @export
rcsbasis <- function(x, knots, derivative = FALSE) {
  if (!is.numeric(x)) stop("`x` must be numeric", call. = FALSE)
  if (!is.numeric(knots) || length(knots) < 2L || any(!is.finite(knots)) || anyDuplicated(knots)) {
    stop("`knots` must hold at least two distinct finite values (the boundary knots)", call. = FALSE)
  }
  rcs_basis_cpp(as.numeric(x), sort(as.numeric(knots)), isTRUE(derivative))
}

#' Default knot placement for a Royston-Parmar spline
#'
#' Places boundary knots at the min/max of `x` and interior knots at
#' equally spaced centiles of `x` restricted to event times, matching the
#' default behaviour of `rstpm2::stpm2()` and `flexsurv::flexsurvspline()`.
#'
#' When event times are heavily concentrated near one end of the range
#' (e.g. an accelerating hazard piling events up before administrative
#' censoring), a centile can coincide with a boundary knot. An interior
#' knot placed exactly at a boundary knot makes that spline basis column
#' identically zero (see `rcsbasis()`), which is a silent rank deficiency:
#' such knots are dropped, reducing the realised spline df below `df` (with
#' a warning) rather than fitting a structurally unidentified model.
#'
#' The knots are identical to those chosen by `rstpm2::stpm2(..., df = df)`
#' for the baseline (centiles of log event time, computed with
#' `stats::quantile()` type 7).
#'
#' @param x numeric vector, typically log(event time) for uncensored observations.
#' @param df degrees of freedom of the spline (number of interior knots + 1).
#' @return a sorted numeric vector of knots, boundary knots first and last,
#'   to pass to [rcsbasis()] or to the `knots` argument of [rpsurv()].
#' @seealso [rcsbasis()].
#' @examples
#' lt <- log(c(5, 10, 20, 40, 80, 160))
#' defaultknots(lt, df = 3)
#' @export
defaultknots <- function(x, df) {
  if (length(df) != 1L || !is.finite(df) || df < 1 || df != round(df)) {
    stop("df must be a single whole number >= 1", call. = FALSE)
  }
  df <- as.integer(df)
  nk <- df - 1L
  if (nk == 0L) {
    probs <- numeric(0)
  } else {
    probs <- seq_len(nk) / (nk + 1L)
  }
  kmin <- min(x); kmax <- max(x)
  interior <- if (length(probs)) stats::quantile(x, probs = probs, names = FALSE) else numeric(0)

  eps <- sqrt(.Machine$double.eps) * max(1, kmax - kmin)
  degenerate <- interior < kmin + eps | interior > kmax - eps
  if (any(degenerate)) {
    warning(sprintf(
      "%d interior knot(s) coincide with a boundary knot (event times are concentrated near %s); dropping them, effective spline df is %d instead of %d",
      sum(degenerate), if (mean(interior[degenerate] > kmax - eps) > 0.5) "the maximum" else "the minimum",
      df - sum(degenerate), df
    ), call. = FALSE)
    interior <- interior[!degenerate]
  }
  if (anyDuplicated(interior)) {
    n_dup <- sum(duplicated(interior))
    warning(sprintf(
      "%d interior knot(s) coincide with another knot because of tied values; dropping them, effective spline df is %d instead of %d",
      n_dup, df - sum(degenerate) - n_dup, df
    ), call. = FALSE)
    interior <- unique(interior)
  }

  sort(c(kmin, interior, kmax))
}
