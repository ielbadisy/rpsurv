#' Residuals for a Royston-Parmar model
#'
#' @param object a fitted `"rpsurv"` object.
#' @param type `"coxsnell"` (Cox-Snell residuals, `r_i = H(stop_i) - H(entry_i)`),
#'   `"martingale"` (`status_i - r_i`), or `"deviance"` (the usual signed
#'   transform of the martingale residual, roughly symmetric around 0 for a
#'   well-fitting model).
#' @param ... unused.
#' @return numeric vector of residuals, one per observation used to fit the model.
#' @importFrom stats residuals
#' @export
residuals.rpsurv <- function(object, type = c("coxsnell", "martingale", "deviance"), ...) {
  type <- match.arg(type)
  eta <- as.numeric(object$X %*% object$coefficients)
  H <- rp_cumhaz(eta, object$scale)
  if (any(object$has_entry > 0)) {
    eta_e <- as.numeric(object$Xentry %*% object$coefficients)
    H <- H - object$has_entry * rp_cumhaz(eta_e, object$scale)
  }
  r <- H
  n_bad <- sum(!is.finite(r) | r < 0)
  if (n_bad) {
    warning(sprintf(
      "%d Cox-Snell residual(s) are negative or not finite: the fitted cumulative hazard is not increasing over those subjects' follow-up, so the model is invalid there",
      n_bad
    ), call. = FALSE)
  }

  if (type == "coxsnell") return(r)

  m <- object$status - r
  if (type == "martingale") return(m)

  dev <- sign(m) * sqrt(pmax(-2 * (m + object$status * log(pmax(r, 1e-12))), 0))
  dev
}

#' Cox-Snell residual diagnostic plot
#'
#' Plots the Nelson-Aalen cumulative hazard of the Cox-Snell residuals
#' against the residuals themselves. Under a correctly specified model the
#' residuals are a censored sample from a unit exponential distribution,
#' so the points should lie close to the y = x reference line.
#' `plot(fit, type = "coxsnell")` draws the same plot.
#'
#' For counting-process data (`Surv(start, stop, status)`), each row's
#' residual is `H(stop) - H(start)`. Given survival to `start`, it is unit
#' exponential and censored at its observed value when the row has no
#' event, so the same check applies row by row. Tied residuals are handled
#' by the Nelson-Aalen estimator.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param ... further arguments passed to [graphics::plot()]. Both axes
#'   default to the same range starting at 0, so the reference line is the
#'   diagonal of the plot.
#' @return Called for its side effect of drawing the diagnostic plot.
#'   Invisibly returns a `data.frame` with columns `time` (the sorted
#'   distinct Cox-Snell residuals) and `cumhaz` (their Nelson-Aalen
#'   cumulative hazard), the coordinates of the plotted step function.
#' @export
coxsnellplot <- function(object, ...) {
  r <- residuals(object, type = "coxsnell")
  ok <- is.finite(r)
  fit <- survival::survfit(survival::Surv(r[ok], object$status[ok]) ~ 1, ctype = 1)
  H <- fit$cumhaz
  lim <- c(0, max(c(fit$time, H), na.rm = TRUE))
  dots <- list(...)
  def <- list(xlim = lim, ylim = lim)
  rp_empty_plot(lim, lim, "Cox-Snell residual", "Nelson-Aalen cumulative hazard of residuals",
                c(def[setdiff(names(def), names(dots))], dots))
  graphics::abline(0, 1, col = "red", lty = 2)
  graphics::lines(c(0, fit$time), c(0, H), type = "s", lwd = 2)
  invisible(data.frame(time = fit$time, cumhaz = H))
}

#' Compare the fitted survival curve against the Kaplan-Meier estimate
#'
#' A calibration / goodness-of-fit diagnostic: overlays the model-based
#' survival curve on the nonparametric Kaplan-Meier estimate, overall or by
#' groups of a covariate. Close agreement supports the chosen spline df and
#' scale; systematic divergence suggests more baseline df, a different
#' `scale`, or a `tve` term is needed. `plot(fit, type = "km")` draws the
#' same plot.
#'
#' The Kaplan-Meier curve estimates the marginal survival of the group, so
#' the model curve is the matching population-averaged (direct-adjusted)
#' estimand, not the survival curve at the mean covariate profile (which is
#' a different, generally higher or lower, curve when covariates vary):
#'
#' * Right-censored data: the mean over the group's subjects of their
#'   predicted survival curves, `mean_i S(t | x_i)`. Kaplan-Meier estimates
#'   the same curve when censoring does not depend on the covariates; when
#'   it does, the two can differ slightly even for a correct model.
#' * Counting-process data (`Surv(start, stop, status)`, left truncation or
#'   time-varying covariates): the Kaplan-Meier curve is then a product of
#'   conditional survival over the rows at risk, so the model curve is
#'   `exp(-A(t))`, where `A` cumulates the average model hazard of the rows
#'   at risk at each time (the model counterpart of the Nelson-Aalen
#'   estimator, computed on a fine time grid). A mean of `S(t | x_i)` from
#'   time 0 would not be comparable, because subjects enter the sample only
#'   after surviving to their entry time.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param by optional grouping: the name of a variable of the model
#'   formula, or a vector with one value per observation used in the fit
#'   (rows with missing values dropped), for example `cut(age, 3)`. At most
#'   10 groups. Each group gets its Kaplan-Meier curve and its own model
#'   curve, averaged over the subjects of that group.
#' @param col colours. Without `by`, two colours for the Kaplan-Meier and the
#'   model curve (default black and red); with `by`, one colour per group
#'   (both curves of a group share it and differ by line type).
#' @param n_grid number of time points at which the model curve is
#'   evaluated, from 0 to the last follow-up time.
#' @param legend `TRUE` (default; the legend goes in the plot corner the
#'   curves cross least), `FALSE` for none, or a position keyword accepted
#'   by [graphics::legend()].
#' @param ... further arguments passed to [graphics::plot()].
#' @return Called for its side effect of drawing the comparison plot.
#'   Invisibly returns a `data.frame` with columns `group`, `curve`
#'   (`"km"` or `"model"`), `time` and `surv`, the plotted coordinates.
#' @examples
#' \donttest{
#' fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype, data = veteran)
#' kmcompareplot(fit, by = "celltype")
#' }
#' @export
kmcompareplot <- function(object, by = NULL, col = NULL, n_grid = 200,
                            legend = TRUE, ...) {
  n <- length(object$time)
  by_name <- NULL
  if (is.null(by)) {
    g <- factor(rep("All", n))
  } else {
    if (length(by) == 1L && is.character(by) && n != 1L) {
      by_name <- by
      avail <- union(names(object$cov_vars), names(object$cov_data))
      v <- if (!is.null(object$cov_vars) && by %in% names(object$cov_vars)) {
        object$cov_vars[[by]]
      } else if (by %in% names(object$cov_data)) {
        object$cov_data[[by]]
      } else if (!length(avail)) {
        stop("`by` names a variable, but the model has no covariates; pass a grouping vector instead",
             call. = FALSE)
      } else {
        stop("`by` must name a variable of the model formula (", paste(avail, collapse = ", "),
             ") or be a grouping vector with one value per fitted observation", call. = FALSE)
      }
    } else if (length(by) == n) {
      v <- by
    } else {
      stop(sprintf("`by` must be a variable name or a grouping vector of length %d (one value per fitted observation)", n),
           call. = FALSE)
    }
    k <- length(unique(v[!is.na(v)]))
    if (k > 10L) {
      stop(sprintf("`by` has %d distinct values; group it into at most 10 (for example `by = cut(x, 3)` with a grouping vector)",
                   k), call. = FALSE)
    }
    g <- if (is.factor(v)) droplevels(v) else factor(v)
  }
  groups <- split(seq_len(n), g)
  groups <- groups[lengths(groups) > 0L]
  k <- length(groups)
  if (is.null(col)) col <- if (k == 1L) c("black", "red") else grDevices::hcl.colors(k, "Dark 3")
  col_km <- if (k == 1L) col[1L] else rep_len(col, k)
  col_md <- if (k == 1L) col[min(2L, length(col))] else col_km

  counting <- any(object$has_entry > 0)
  tmax <- max(object$time)
  grid <- seq(0, tmax, length.out = n_grid + 1L)

  rp_empty_plot(c(0, tmax), c(0, 1), "Time", "Survival probability", list(...))
  out <- vector("list", k)
  for (i in seq_len(k)) {
    idx <- groups[[i]]
    sf <- if (counting) {
      survival::survfit(survival::Surv(object$entry[idx], object$time[idx], object$status[idx]) ~ 1)
    } else {
      survival::survfit(survival::Surv(object$time[idx], object$status[idx]) ~ 1)
    }
    km_t <- c(0, sf$time)
    km_s <- c(1, sf$surv)
    graphics::lines(km_t, km_s, type = "s", col = col_km[i], lty = 1, lwd = 1.5)

    s_model <- rp_average_survival(object, idx, grid, counting)
    graphics::lines(grid, s_model, col = col_md[i], lty = 2, lwd = 2)
    nm <- names(groups)[i]
    out[[i]] <- rbind(data.frame(group = nm, curve = "km", time = km_t, surv = km_s),
                      data.frame(group = nm, curve = "model", time = grid, surv = s_model))
  }

  if (isTRUE(legend)) {
    pts <- do.call(rbind, out)
    legend <- rp_legend_pos(pts$time, pts$surv, if (k == 1L) 2L else k + 2L)
  }
  if (!isFALSE(legend)) {
    if (k == 1L) {
      graphics::legend(legend, legend = c("Kaplan-Meier", "rpsurv (population-averaged)"),
                       col = c(col_km, col_md), lty = c(1, 2), lwd = c(1.5, 2), bty = "n",
                       inset = 0.02)
    } else {
      graphics::legend(legend, legend = c(if (is.null(by_name)) names(groups) else paste(by_name, "=", names(groups)),
                                    "Kaplan-Meier", "rpsurv"),
                       col = c(col_km, "black", "black"), lty = c(rep(NA, k), 1, 2),
                       pch = c(rep(15, k), NA, NA), lwd = c(rep(NA, k), 1.5, 2),
                       pt.cex = 1.5, bty = "n", inset = 0.02)
    }
  }
  invisible(do.call(rbind, out))
}

# Model survival averaged over the fitted rows `rows`, on `grid` (starting
# at 0). Right-censored data: mean of the individual survival curves.
# Counting-process data: exp(-A(t)), A cumulating over each grid interval
# the summed cumulative-hazard increments of the rows at risk divided by the
# mean number at risk over the interval. Computed in chunks of rows.
#' @keywords internal
rp_average_survival <- function(object, rows, grid, counting) {
  tg <- grid[-1L]
  G <- length(tg)
  b <- object$coefficients
  lt <- log(tg)
  if (counting) {
    H_exit <- rp_cumhaz(as.numeric(object$X %*% b), object$scale)
    H_entry <- object$has_entry * rp_cumhaz(as.numeric(object$Xentry %*% b), object$scale)
    sum_dH <- numeric(G)
    sum_len <- numeric(G)
  } else {
    sum_S <- numeric(G)
  }
  chunk <- max(1L, floor(2e5 / G))
  for (s in seq(1L, length(rows), by = chunk)) {
    r <- rows[s:min(length(rows), s + chunk - 1L)]
    m <- length(r)
    cov <- object$cov_data[rep(r, each = G), , drop = FALSE]
    X <- rp_design(rep(lt, m), cov, object$knots, object$tve, object$tve_knots,
                   object$nle, object$nle_knots)$X
    H <- matrix(rp_cumhaz(as.numeric(X %*% b), object$scale), nrow = m, byrow = TRUE)
    if (!counting) {
      sum_S <- sum_S + colSums(exp(-H))
      next
    }
    H <- cbind(0, H) # H at every grid point, including time 0
    lo_t <- matrix(grid[-(G + 1L)], m, G, byrow = TRUE)
    hi_t <- matrix(tg, m, G, byrow = TRUE)
    en <- object$entry[r]
    ex <- object$time[r]
    H_hi <- ifelse(hi_t <= ex, H[, -1L, drop = FALSE], H_exit[r])
    H_lo <- ifelse(lo_t >= en, H[, -(G + 1L), drop = FALSE], H_entry[r])
    len <- pmin(hi_t, ex) - pmax(lo_t, en) # matrices: keep the dim attribute
    len[len < 0] <- 0
    dH <- ifelse(len > 0, H_hi - H_lo, 0)
    sum_dH <- sum_dH + colSums(dH)
    sum_len <- sum_len + colSums(len)
  }
  if (!counting) return(c(1, sum_S / length(rows)))
  inc <- ifelse(sum_len > 0, sum_dH / (sum_len / diff(grid)), 0)
  c(1, exp(-cumsum(inc)))
}
