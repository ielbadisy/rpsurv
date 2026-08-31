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
  inv_link <- rp_inv_link(object$scale)

  eta <- as.numeric(object$X %*% object$coefficients)
  H <- -log(inv_link(eta))
  if (any(object$has_entry > 0)) {
    eta_e <- as.numeric(object$Xentry %*% object$coefficients)
    H <- H - object$has_entry * (-log(inv_link(eta_e)))
  }
  r <- H

  if (type == "coxsnell") return(r)

  m <- object$status - r
  if (type == "martingale") return(m)

  dev <- sign(m) * sqrt(pmax(-2 * (m + object$status * log(pmax(r, 1e-12))), 0))
  dev
}

#' Cox-Snell residual diagnostic plot
#'
#' Plots the Nelson-Aalen cumulative hazard of the Cox-Snell residuals
#' against the residuals themselves. A well-fitting model gives points
#' scattered around the y = x line.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param ... further arguments passed to [graphics::plot()].
#' @return Called for its side effect of drawing the diagnostic plot.
#'   Invisibly returns a `data.frame` with columns `time` (the sorted
#'   Cox-Snell residuals) and `cumhaz` (their Nelson-Aalen cumulative
#'   hazard), the coordinates of the plotted step function.
#' @export
coxsnell_plot <- function(object, ...) {
  r <- residuals(object, type = "coxsnell")
  fit <- survival::survfit(survival::Surv(r, object$status) ~ 1)
  H <- -log(fit$surv)
  graphics::plot(fit$time, H, xlab = "Cox-Snell residual", ylab = "Cumulative hazard of residuals",
                 type = "s", ...)
  graphics::abline(0, 1, col = "red", lty = 2)
  invisible(data.frame(time = fit$time, cumhaz = H))
}

#' Compare the fitted survival curve against the Kaplan-Meier estimate
#'
#' A calibration / goodness-of-fit diagnostic: overlays the model-predicted
#' survival curve on the nonparametric Kaplan-Meier estimate, optionally by
#' strata of a categorical covariate. Close agreement supports the chosen
#' spline df and scale; systematic divergence suggests more baseline df,
#' a different `scale`, or a `tve` term is needed.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param by optional name of a covariate in the fitted data to stratify by
#'   (each level gets its own KM curve and predicted curve, at that level's
#'   mean of the other covariates).
#' @param col colours for KM (solid step) vs model (dashed) curves.
#' @param ... further arguments passed to [graphics::plot()].
#' @return No return value, called for its side effect of drawing the
#'   Kaplan-Meier versus fitted-survival comparison plot.
#' @export
km_compare_plot <- function(object, by = NULL, col = c("black", "red"), ...) {
  d <- data.frame(time = object$time, status = object$status, entry = object$entry, object$cov_data)

  groups <- if (is.null(by)) list(seq_len(nrow(d))) else split(seq_len(nrow(d)), d[[by]])
  group_names <- if (is.null(by)) "overall" else names(groups)

  graphics::plot(range(d$time[d$time > 0]), c(0, 1), type = "n",
                 xlab = "Time", ylab = "Survival probability", ...)

  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    sf <- if (any(d$entry[idx] > 0)) {
      survival::survfit(survival::Surv(entry, time, status) ~ 1, data = d[idx, , drop = FALSE])
    } else {
      survival::survfit(survival::Surv(time, status) ~ 1, data = d[idx, , drop = FALSE])
    }
    graphics::lines(sf$time, sf$surv, type = "s", col = col[1], lty = 1)

    newdata <- as.data.frame(lapply(object$cov_data[idx, , drop = FALSE], function(v) mean(as.numeric(v))))
    if (!is.null(by)) newdata[[by]] <- d[[by]][idx[1]]
    pred <- predict(object, newdata = newdata, type = "survival")
    graphics::lines(pred$time, pred$est, col = col[2], lty = 2, lwd = 2)
  }

  graphics::legend("bottomleft", legend = c("Kaplan-Meier", "rpsurv"), col = col, lty = c(1, 2), bty = "n")
  invisible(NULL)
}
