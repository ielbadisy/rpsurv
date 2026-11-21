#' Plot a Royston-Parmar model
#'
#' Draws predicted survival, hazard or cumulative hazard curves for
#' covariate profiles, or one of the two goodness-of-fit diagnostics:
#' `type = "coxsnell"` calls [coxsnellplot()] and `type = "km"` calls
#' [kmcompareplot()].
#'
#' @param x a fitted `"rpsurv"` object.
#' @param newdata covariate profile(s) to plot, one row per curve, on the
#'   original variables of the model formula (see [predict.rpsurv()]).
#'   Defaults to a single profile with each numeric variable at its mean
#'   and each factor, character or logical variable at its most common
#'   value in the fitted data.
#' @param type one of `"survival"`, `"hazard"`, `"cumhaz"` (predicted
#'   curves), `"coxsnell"` (Cox-Snell residual plot, see [coxsnellplot()])
#'   or `"km"` (Kaplan-Meier vs model, see [kmcompareplot()]). The
#'   diagnostics use the fitted data, so `newdata`, `times` and `ci` do not
#'   apply to them.
#' @param times times at which to evaluate the curve(s). Defaults to 100
#'   points from the first event time to the last follow-up time, the range
#'   over which the model is estimated.
#' @param ci logical; add a 95% delta-method confidence band (pointwise).
#' @param col colour(s), recycled over rows of `newdata`. By default black,
#'   blue, vermillion, then the other Okabe-Ito colours (colour-blind safe).
#'   For `type = "coxsnell"`, the colour of the residual curve.
#' @param legend `TRUE`, `FALSE`, or a position keyword accepted by
#'   [graphics::legend()] for the legend identifying the curves, labelled
#'   by the `newdata` columns that differ between rows. By default a legend
#'   is drawn when there are several curves, in the plot corner that the
#'   curves cross least. For `type = "km"` the legend is drawn by default.
#'   For predicted curves, a character vector with one label per row of
#'   `newdata` (that is not a single position keyword) is taken as the
#'   labels, shorthand for `labels =`.
#' @param by for `type = "km"` only: optional grouping, the name of a
#'   variable of the model formula or a vector with one value per
#'   observation used in the fit (see [kmcompareplot()]).
#' @param labels optional character vector, one legend label per row of
#'   `newdata` (predicted curves only), replacing the default labels built
#'   from the `newdata` columns, e.g. `labels = c("Remission", "Relapse")`.
#'   Giving `labels` draws the legend unless `legend = FALSE`. For
#'   `type = "km"`, the two curve labels of [kmcompareplot()].
#' @param ... further arguments passed to [graphics::plot()], such as
#'   `ylim` or `main` (and `n_grid`, `minn`, `levels`, `grouplabels` for
#'   `type = "km"`, `ref.col` for `type = "coxsnell"`). No title is drawn
#'   unless `main` is given.
#' @return Called for its side effect of drawing the plot. For the
#'   diagnostics, invisibly returns the value of [coxsnellplot()] or
#'   [kmcompareplot()]. For predicted curves, invisibly returns the
#'   `data.frame` of predictions produced by [predict.rpsurv()], with columns `id` (row of `newdata`), `time`,
#'   `est`, and, when `ci = TRUE`, `lower` and `upper` for the 95%
#'   confidence band. These are exactly the plotted values.
#' @examples
#' \donttest{
#' fit <- rpsurv(survival::Surv(time, status) ~ trt + celltype, data = veteran)
#' plot(fit, newdata = data.frame(trt = 1, celltype = c("squamous", "adeno")))
#' plot(fit, type = "coxsnell")
#' plot(fit, type = "km", by = "celltype")
#' }
#' @export
plot.rpsurv <- function(x, newdata = NULL,
                        type = c("survival", "hazard", "cumhaz", "coxsnell", "km"),
                        times = NULL, ci = TRUE, col = NULL, legend = NULL, by = NULL,
                        labels = NULL, ...) {
  type <- match.arg(type)
  if (type %in% c("coxsnell", "km")) {
    unused <- c(newdata = !is.null(newdata), times = !is.null(times), ci = !missing(ci))
    if (type == "coxsnell") unused <- c(unused, legend = !is.null(legend), labels = !is.null(labels))
    if (any(unused)) {
      stop("`", paste(names(unused)[unused], collapse = "`, `"), "` ",
           if (sum(unused) > 1L) "do" else "does", " not apply to type = \"", type,
           "\": the diagnostic uses the fitted data.", call. = FALSE)
    }
    if (type == "coxsnell") {
      if (!is.null(by)) stop("`by` applies to type = \"km\" only.", call. = FALSE)
      return(if (is.null(col)) coxsnellplot(x, ...) else coxsnellplot(x, col = col, ...))
    }
    km <- list(x, by = by, col = col, legend = if (is.null(legend)) TRUE else legend)
    return(do.call(kmcompareplot, c(km, if (!is.null(labels)) list(labels = labels), list(...))))
  }
  if (!is.null(by)) stop("`by` applies to type = \"km\" only.", call. = FALSE)
  if (is.null(newdata)) newdata <- rp_default_profile(x)
  newdata <- as.data.frame(newdata)
  col <- rep_len(if (is.null(col)) rp_palette(max(1L, nrow(newdata))) else col, max(1L, nrow(newdata)))

  pred <- predict(x, newdata = newdata, times = times, type = type, se.fit = ci)
  ylab <- switch(type, survival = "Survival probability", hazard = "Hazard", cumhaz = "Cumulative hazard")
  rp_draw_curves(pred, newdata, ylab, type == "survival", col, legend, list(...), labels)
}

# Curves (and 95% bands when `pred` has them) of a long prediction data frame, one per row of `newdata`;
# `prob` fixes the y range to [0, 1]. Returns `pred` invisibly.
#' @keywords internal
rp_draw_curves <- function(pred, newdata, ylab, prob, col, legend, dots, labels = NULL) {
  n_curves <- max(1L, nrow(newdata))
  lg <- rp_legend_spec(legend, labels, newdata); legend <- lg$legend; labels <- lg$labels
  has_band <- !is.null(pred$lower)
  yr <- if (prob) c(0, 1) else range(c(0, pred$est, if (has_band) c(pred$lower, pred$upper)), finite = TRUE)
  rp_empty_plot(range(pred$time), yr, "Time", ylab, dots)
  for (i in seq_len(n_curves)) {
    d <- pred[pred$id == i, ]
    if (has_band) {
      ok <- is.finite(d$lower) & is.finite(d$upper)
      graphics::polygon(c(d$time[ok], rev(d$time[ok])), c(d$lower[ok], rev(d$upper[ok])),
                         col = grDevices::adjustcolor(col[i], alpha.f = 0.15), border = NA)
    }
    graphics::lines(d$time, d$est, col = col[i], lwd = 2)
  }

  if (isTRUE(legend)) {
    legend <- rp_legend_pos(rep(pred$time, if (has_band) 3L else 1L),
                            c(pred$est, if (has_band) c(pred$lower, pred$upper)), n_curves)
  }
  if (!isFALSE(legend)) {
    graphics::legend(legend, legend = labels, col = col, lwd = 2,
                     bty = "n", inset = 0.02)
  }
  invisible(pred)
}

# Default colours: black, blue, vermillion, then the rest of Okabe-Ito (colour-blind safe), recycled.
#' @keywords internal
rp_palette <- function(k) {
  rep_len(c("black", "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442",
            "grey50", "grey75"), k)
}

# Legend position (TRUE = automatic, FALSE = none) and labels. A character `legend` that is not a single
# position keyword is the labels; custom labels switch the legend on.
#' @keywords internal
rp_legend_spec <- function(legend, labels, newdata) {
  n <- max(1L, nrow(newdata))
  pos <- c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")
  if (is.character(legend) && !(length(legend) == 1L && legend %in% pos)) {
    labels <- legend; legend <- TRUE
  }
  if (!is.null(labels) && length(labels) != n) {
    stop("`labels` must have one entry per row of `newdata` (", n, ").", call. = FALSE)
  }
  if (is.null(legend)) legend <- n > 1L || !is.null(labels)
  list(legend = legend, labels = if (is.null(labels)) rp_profile_labels(newdata) else as.character(labels))
}

# One covariate profile: numeric variables at their mean, factor, character
# and logical variables at their most common value in the fitted data.
#' @keywords internal
rp_default_profile <- function(object) {
  if (!ncol(object$cov_data)) return(data.frame(row.names = 1L))
  vars <- object$cov_vars
  if (is.null(vars)) {
    # fits without the original variables: mean of the design columns
    return(as.data.frame(lapply(object$cov_data, function(v) mean(as.numeric(v))),
                         optional = TRUE))
  }
  prof <- lapply(vars, function(v) {
    if (is.numeric(v) && !is.factor(v)) return(mean(v, na.rm = TRUE))
    tab <- table(v)
    m <- names(tab)[which.max(tab)]
    if (is.factor(v)) factor(m, levels = levels(v)) else if (is.logical(v)) as.logical(m) else m
  })
  as.data.frame(prof, optional = TRUE, stringsAsFactors = FALSE)
}

# Legend labels for the rows of `newdata`, from the columns that differ
# between rows ("trt = 1, celltype = adeno").
#' @keywords internal
rp_profile_labels <- function(newdata) {
  n <- max(1L, nrow(newdata))
  if (!ncol(newdata)) return(as.character(seq_len(n)))
  fmt <- function(v) if (is.numeric(v)) format(signif(v, 3), trim = TRUE) else as.character(v)
  vary <- vapply(newdata, function(v) length(unique(v)) > 1L, logical(1))
  if (!any(vary)) vary[] <- TRUE
  cols <- lapply(names(newdata)[vary], function(nm) paste(nm, "=", fmt(newdata[[nm]])))
  do.call(paste, c(cols, sep = ", "))
}

# Empty plot frame whose default axis labels can be overridden through `...`
# (passing `xlab` or `ylab` there would otherwise match twice).
#' @keywords internal
rp_empty_plot <- function(xr, yr, xlab, ylab, dots) {
  def <- list(xlab = xlab, ylab = ylab)
  dots <- c(def[setdiff(names(def), names(dots))], dots)
  do.call(graphics::plot, c(list(xr, yr, type = "n"), dots))
}

# Corner of the plot region crossed by the fewest curve points, for a legend
# of `n_entries` lines: each candidate box spans 45% of the width and a
# height growing with the number of entries.
#' @keywords internal
rp_legend_pos <- function(x, y, n_entries) {
  ok <- is.finite(x) & is.finite(y)
  usr <- graphics::par("usr")
  if (graphics::par("ylog")) y <- log10(y)
  if (graphics::par("xlog")) x <- log10(x)
  u <- (x[ok] - usr[1L]) / (usr[2L] - usr[1L])
  v <- (y[ok] - usr[3L]) / (usr[4L] - usr[3L])
  h <- min(0.9, 0.08 + 0.065 * n_entries)
  w <- 0.45
  pos <- c("topright", "bottomleft", "topleft", "bottomright")
  hits <- c(sum(u > 1 - w & v > 1 - h), sum(u < w & v < h),
            sum(u < w & v > 1 - h), sum(u > 1 - w & v < h))
  pos[which.min(hits)]
}
