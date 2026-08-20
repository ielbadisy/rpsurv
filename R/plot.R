#' Plot predicted curves from a Royston-Parmar model
#'
#' @param x a fitted `"rpsurv"` object.
#' @param newdata covariate profile(s) to plot, one row per curve. Defaults
#'   to a single profile at the mean of each covariate.
#' @param type one of `"survival"`, `"hazard"`, `"cumhaz"`.
#' @param times times at which to evaluate the curve(s).
#' @param ci logical; add a 95% confidence band (ignored for `type = "hazard"`).
#' @param col colour(s), recycled over rows of `newdata`.
#' @param ... further arguments passed to [graphics::plot()].
#' @export
plot.rpsurv <- function(x, newdata = NULL, type = c("survival", "hazard", "cumhaz"),
                         times = NULL, ci = TRUE, col = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (ncol(x$cov_data)) {
      newdata <- as.data.frame(lapply(x$cov_data, function(v) mean(as.numeric(v))))
    } else {
      newdata <- data.frame(row.names = 1L)
    }
  }
  n_curves <- max(1L, nrow(newdata))
  if (is.null(col)) col <- seq_len(n_curves)

  pred <- predict(x, newdata = newdata, times = times, type = type, se.fit = ci && type != "hazard")

  ylab <- switch(type, survival = "Survival probability", hazard = "Hazard", cumhaz = "Cumulative hazard")
  graphics::plot(range(pred$time), range(pred$est), type = "n", xlab = "Time", ylab = ylab, ...)
  for (i in seq_len(n_curves)) {
    d <- pred[pred$id == i, ]
    if (ci && type != "hazard" && !is.null(d$lower)) {
      graphics::polygon(c(d$time, rev(d$time)), c(d$lower, rev(d$upper)),
                         col = grDevices::adjustcolor(col[i], alpha.f = 0.15), border = NA)
    }
    graphics::lines(d$time, d$est, col = col[i], lwd = 2)
  }
  invisible(pred)
}
