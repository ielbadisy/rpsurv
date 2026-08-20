#' Residuals for a Royston-Parmar model
#'
#' @param object a fitted `"rpsurv"` object.
#' @param type `"coxsnell"` (Cox-Snell residuals, `r_i = H(t_i | x_i)`) or
#'   `"martingale"` (`d_i - r_i`).
#' @param ... unused.
#' @return numeric vector of residuals, one per observation used to fit the model.
#' @importFrom stats residuals
#' @export
residuals.rpsurv <- function(object, type = c("coxsnell", "martingale"), ...) {
  type <- match.arg(type)
  inv_link <- rp_inv_link(object$scale)
  eta <- as.numeric(object$X %*% object$coefficients)
  S <- inv_link(eta)
  r <- -log(S)
  if (type == "coxsnell") return(r)
  object$status - r
}

#' Cox-Snell residual diagnostic plot
#'
#' Plots the Nelson-Aalen cumulative hazard of the Cox-Snell residuals
#' against the residuals themselves. A well-fitting model gives points
#' scattered around the y = x line.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param ... further arguments passed to [graphics::plot()].
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
