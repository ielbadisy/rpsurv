#' Restricted cubic spline basis (Durrleman-Simon / Royston-Parmar parameterisation)
#'
#' Builds the basis used by Royston & Parmar (2002) flexible parametric
#' models: a linear term plus one term per interior knot, each term being
#' a natural (restricted) cubic spline component that is linear beyond the
#' boundary knots.
#'
#' @param x numeric vector (typically log time) at which to evaluate the basis.
#' @param knots numeric vector of ALL knots (boundary knots first and last,
#'   interior knots in between), already sorted.
#' @param derivative if `TRUE`, return d(basis)/dx instead of the basis itself.
#' @return a matrix with `length(knots) - 1` columns.
#' @keywords internal
rcs_basis <- function(x, knots, derivative = FALSE) {
  knots <- sort(knots)
  kmin <- knots[1L]
  kmax <- knots[length(knots)]
  interior <- knots[-c(1L, length(knots))]
  nk <- length(interior)

  pos <- function(u) {
    p <- u
    p[u < 0] <- 0
    p
  }

  if (!derivative) {
    out <- matrix(0, length(x), nk + 1L)
    out[, 1L] <- x
    for (j in seq_len(nk)) {
      kj <- interior[j]
      lambda_j <- (kmax - kj) / (kmax - kmin)
      out[, j + 1L] <- pos(x - kj)^3 -
        lambda_j * pos(x - kmin)^3 -
        (1 - lambda_j) * pos(x - kmax)^3
    }
  } else {
    out <- matrix(0, length(x), nk + 1L)
    out[, 1L] <- 1
    for (j in seq_len(nk)) {
      kj <- interior[j]
      lambda_j <- (kmax - kj) / (kmax - kmin)
      out[, j + 1L] <- 3 * pos(x - kj)^2 -
        lambda_j * 3 * pos(x - kmin)^2 -
        (1 - lambda_j) * 3 * pos(x - kmax)^2
    }
  }
  out
}

#' Default knot placement for a Royston-Parmar spline
#'
#' Places boundary knots at the min/max of `x` and interior knots at
#' equally spaced centiles of `x` restricted to event times, matching the
#' default behaviour of `rstpm2::stpm2()` / `flexsurv::flexsurvspline()`.
#'
#' @param x numeric vector, typically log(event time) for uncensored observations.
#' @param df degrees of freedom of the spline (number of interior knots + 1).
#' @keywords internal
default_knots <- function(x, df) {
  if (df < 1L) stop("df must be >= 1", call. = FALSE)
  nk <- df - 1L
  if (nk == 0L) {
    probs <- numeric(0)
  } else {
    probs <- seq_len(nk) / (nk + 1L)
  }
  interior <- if (length(probs)) stats::quantile(x, probs = probs, names = FALSE) else numeric(0)
  sort(c(min(x), interior, max(x)))
}
