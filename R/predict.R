rp_inv_link <- function(scale) {
  switch(scale,
    hazard = function(eta) exp(-exp(eta)),
    odds = function(eta) 1 / (1 + exp(eta)),
    normal = function(eta) stats::pnorm(-eta)
  )
}

#' Predict from a fitted Royston-Parmar model
#'
#' @param object a fitted `"rpsurv"` object.
#' @param newdata a data frame of covariate values, one row per subject.
#'   Defaults to the covariate values used to fit the model.
#' @param times numeric vector of times at which to predict. Defaults to
#'   a grid over the observed follow-up.
#' @param type one of `"survival"`, `"hazard"`, `"cumhaz"` or `"link"`.
#' @param se.fit logical; if `TRUE`, adds `lower`/`upper` 95% confidence
#'   limits (via the delta method on the linear predictor scale). Only
#'   supported for `type %in% c("survival", "cumhaz", "link")`.
#' @param ... unused.
#' @return a long-format data frame with columns `id`, `time`, `est`,
#'   and optionally `lower`/`upper`.
#' @importFrom stats predict
#' @export
predict.rpsurv <- function(object, newdata = NULL, times = NULL,
                            type = c("survival", "hazard", "cumhaz", "link"),
                            se.fit = FALSE, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) newdata <- object$cov_data
  newdata <- as.data.frame(newdata)
  if (ncol(object$cov_data) && !all(colnames(object$cov_data) %in% colnames(newdata))) {
    stop("`newdata` is missing covariate(s): ",
         paste(setdiff(colnames(object$cov_data), colnames(newdata)), collapse = ", "), call. = FALSE)
  }

  if (is.null(times)) {
    times <- seq(min(object$time[object$status == 1]), max(object$time), length.out = 100)
  }
  log_times <- log(times)

  inv_link <- rp_inv_link(object$scale)
  b <- object$coefficients
  V <- object$vcov

  results <- vector("list", nrow(newdata))
  for (i in seq_len(nrow(newdata))) {
    cov_i <- newdata[i, colnames(object$cov_data), drop = FALSE]
    cov_rep <- cov_i[rep(1L, length(log_times)), , drop = FALSE]
    design <- rp_design(log_times, cov_rep, object$knots, object$tve, object$tve_knots)
    X <- design$X
    dX <- design$dX
    eta <- as.numeric(X %*% b)
    deta <- as.numeric(dX %*% b)

    est <- switch(type,
      link = eta,
      survival = inv_link(eta),
      cumhaz = -log(inv_link(eta)),
      hazard = {
        S <- inv_link(eta)
        h <- switch(object$scale,
          hazard = deta / times * exp(eta),
          odds = deta / times * exp(eta) * S,
          normal = deta / times * exp(stats::dnorm(eta, log = TRUE)) / S
        )
        h
      }
    )

    out <- data.frame(id = i, time = times, est = est)

    if (se.fit && type %in% c("survival", "cumhaz", "link")) {
      se_eta <- sqrt(rowSums((X %*% V) * X))
      lo_eta <- eta - stats::qnorm(0.975) * se_eta
      hi_eta <- eta + stats::qnorm(0.975) * se_eta
      if (type == "link") {
        out$lower <- lo_eta
        out$upper <- hi_eta
      } else if (type == "survival") {
        out$lower <- inv_link(hi_eta)
        out$upper <- inv_link(lo_eta)
      } else {
        out$lower <- -log(inv_link(hi_eta))
        out$upper <- -log(inv_link(lo_eta))
      }
    }
    results[[i]] <- out
  }

  do.call(rbind, results)
}
