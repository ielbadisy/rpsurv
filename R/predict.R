rp_inv_link <- function(scale) {
  switch(scale,
    hazard = function(eta) exp(-exp(eta)),
    odds = function(eta) 1 / (1 + exp(eta)),
    normal = function(eta) stats::pnorm(-eta)
  )
}

#' @keywords internal
rp_predict_design <- function(object, cov_i, log_times) {
  cov_rep <- cov_i[rep(1L, length(log_times)), , drop = FALSE]
  design <- rp_design(log_times, cov_rep, object$knots, object$tve, object$tve_knots)
  X <- design$X
  dX <- design$dX
  eta <- as.numeric(X %*% object$coefficients)
  deta <- as.numeric(dX %*% object$coefficients)
  S <- rp_inv_link(object$scale)(eta)
  list(X = X, dX = dX, eta = eta, deta = deta, S = S)
}

#' @keywords internal
rp_hazard_from_design <- function(d, times, scale) {
  switch(scale,
    hazard = d$deta / times * exp(d$eta),
    odds = d$deta / times * exp(d$eta) * d$S,
    normal = d$deta / times * exp(stats::dnorm(d$eta, log = TRUE)) / d$S
  )
}

# Gradient (wrt beta) of log(hazard), one row per time point.
#' @keywords internal
rp_grad_loghaz <- function(d, scale) {
  base <- d$dX / d$deta
  switch(scale,
    hazard = base + d$X,
    odds = base + d$S * d$X,
    normal = base - d$eta * d$X + (stats::dnorm(d$eta) / d$S) * d$X
  )
}

# Gradient (wrt beta) of the survival function, one row per time point.
#' @keywords internal
rp_grad_surv <- function(d, scale) {
  switch(scale,
    hazard = -exp(d$eta) * d$S * d$X,
    odds = -d$S * (1 - d$S) * d$X,
    normal = -stats::dnorm(d$eta) * d$X
  )
}

# Quadratic-form SE, one value per row of `grad` (n_times x p matrix).
#' @keywords internal
rp_delta_se <- function(grad, V) sqrt(rowSums((grad %*% V) * grad))

#' Predict from a fitted Royston-Parmar model
#'
#' @param object a fitted `"rpsurv"` object.
#' @param newdata a data frame of covariate values, one row per subject.
#'   Defaults to the covariate values used to fit the model. For
#'   `type %in% c("hr", "sdiff")`, this is the *numerator*/contrast
#'   covariate profile.
#' @param newdata0 a data frame of covariate values, the *reference*
#'   profile for `type %in% c("hr", "sdiff")`, with the same number of
#'   rows as `newdata`. Required (and only used) for those two types.
#' @param times numeric vector of times at which to predict. Defaults to
#'   a grid over the observed follow-up.
#' @param type one of `"survival"`, `"hazard"`, `"cumhaz"`, `"link"`,
#'   `"hr"` or `"sdiff"`.
#'   \describe{
#'     \item{`"hr"`}{The model-implied instantaneous hazard ratio,
#'       `hazard(newdata) / hazard(newdata0)`, computed as the ratio of
#'       two `type = "hazard"` predictions. This is the correct contrast
#'       under a time-varying effect (`tve`), where the naive
#'       `exp(eta1 - eta0)` is not the instantaneous hazard ratio in
#'       general.}
#'     \item{`"sdiff"`}{The survival difference,
#'       `survival(newdata) - survival(newdata0)`.}
#'   }
#' @param se.fit logical; if `TRUE`, adds `lower`/`upper` 95% confidence
#'   limits via the delta method: on the linear-predictor scale for
#'   `"survival"`/`"cumhaz"`/`"link"`; on the log scale for `"hazard"`/
#'   `"hr"` (so limits stay positive); on the natural scale for
#'   `"sdiff"`.
#' @param ... unused.
#' @return a long-format data frame with columns `id`, `time`, `est`,
#'   and optionally `lower`/`upper`.
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   time = rexp(200, 0.2), status = rbinom(200, 1, 0.7), x1 = rbinom(200, 1, 0.5)
#' )
#' fit <- rpsurv(survival::Surv(time, status) ~ x1, data = dat, tve = "x1")
#' tt <- seq(0.5, 4, length.out = 20)
#'
#' hr <- predict(fit,
#'   newdata = data.frame(x1 = 1), newdata0 = data.frame(x1 = 0),
#'   times = tt, type = "hr", se.fit = TRUE
#' )
#' head(hr)
#'
#' sd <- predict(fit,
#'   newdata = data.frame(x1 = 1), newdata0 = data.frame(x1 = 0),
#'   times = tt, type = "sdiff", se.fit = TRUE
#' )
#' head(sd)
#' }
#' @importFrom stats predict
#' @export
predict.rpsurv <- function(object, newdata = NULL, newdata0 = NULL, times = NULL,
                            type = c("survival", "hazard", "cumhaz", "link", "hr", "sdiff"),
                            se.fit = FALSE, ...) {
  type <- match.arg(type)
  z <- stats::qnorm(0.975)

  if (is.null(times)) {
    times <- seq(min(object$time[object$status == 1]), max(object$time), length.out = 100)
  }
  log_times <- log(times)
  inv_link <- rp_inv_link(object$scale)
  V <- object$vcov

  if (type %in% c("hr", "sdiff")) {
    if (is.null(newdata) || is.null(newdata0)) {
      stop("`type = \"", type, "\"` requires both `newdata` (contrast profile) and ",
           "`newdata0` (reference profile).", call. = FALSE)
    }
    newdata <- as.data.frame(newdata)
    newdata0 <- as.data.frame(newdata0)
    if (nrow(newdata) != nrow(newdata0)) {
      stop("`newdata` and `newdata0` must have the same number of rows.", call. = FALSE)
    }

    results <- vector("list", nrow(newdata))
    for (i in seq_len(nrow(newdata))) {
      cov1 <- newdata[i, colnames(object$cov_data), drop = FALSE]
      cov0 <- newdata0[i, colnames(object$cov_data), drop = FALSE]
      d1 <- rp_predict_design(object, cov1, log_times)
      d0 <- rp_predict_design(object, cov0, log_times)

      if (type == "hr") {
        h1 <- rp_hazard_from_design(d1, times, object$scale)
        h0 <- rp_hazard_from_design(d0, times, object$scale)
        est <- h1 / h0
        out <- data.frame(id = i, time = times, est = est)
        if (se.fit) {
          grad <- rp_grad_loghaz(d1, object$scale) - rp_grad_loghaz(d0, object$scale)
          se_log <- rp_delta_se(grad, V)
          out$lower <- est * exp(-z * se_log)
          out$upper <- est * exp(z * se_log)
        }
      } else {
        est <- d1$S - d0$S
        out <- data.frame(id = i, time = times, est = est)
        if (se.fit) {
          grad <- rp_grad_surv(d1, object$scale) - rp_grad_surv(d0, object$scale)
          se_diff <- rp_delta_se(grad, V)
          out$lower <- est - z * se_diff
          out$upper <- est + z * se_diff
        }
      }
      results[[i]] <- out
    }
    return(do.call(rbind, results))
  }

  if (is.null(newdata)) newdata <- object$cov_data
  newdata <- as.data.frame(newdata)
  if (ncol(object$cov_data) && !all(colnames(object$cov_data) %in% colnames(newdata))) {
    stop("`newdata` is missing covariate(s): ",
         paste(setdiff(colnames(object$cov_data), colnames(newdata)), collapse = ", "), call. = FALSE)
  }

  results <- vector("list", nrow(newdata))
  for (i in seq_len(nrow(newdata))) {
    cov_i <- newdata[i, colnames(object$cov_data), drop = FALSE]
    d <- rp_predict_design(object, cov_i, log_times)
    X <- d$X; eta <- d$eta

    est <- switch(type,
      link = eta,
      survival = d$S,
      cumhaz = -log(d$S),
      hazard = rp_hazard_from_design(d, times, object$scale)
    )

    out <- data.frame(id = i, time = times, est = est)

    if (se.fit && type %in% c("survival", "cumhaz", "link")) {
      se_eta <- rp_delta_se(X, V)
      lo_eta <- eta - z * se_eta
      hi_eta <- eta + z * se_eta
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
    } else if (se.fit && type == "hazard") {
      grad <- rp_grad_loghaz(d, object$scale)
      se_log <- rp_delta_se(grad, V)
      out$lower <- est * exp(-z * se_log)
      out$upper <- est * exp(z * se_log)
    }
    results[[i]] <- out
  }

  do.call(rbind, results)
}

#' Restricted mean survival time difference between two covariate profiles
#'
#' Integrates the survival difference `sdiff = S(t | newdata) -
#' S(t | newdata0)` from 0 to `tau` by the trapezoidal rule, with an
#' optional delta-method standard error and 95% confidence interval.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param newdata a one-row data frame, the contrast covariate profile.
#' @param newdata0 a one-row data frame, the reference covariate profile.
#' @param tau numeric restriction time.
#' @param se.fit logical; if `TRUE`, adds a delta-method standard error
#'   and 95% confidence interval. The delta method is applied directly to
#'   the trapezoidal-weighted sum of the survival-difference gradients
#'   (a linear combination of the model coefficients' asymptotic normal
#'   distribution stays normal), not by resampling.
#' @param n_grid number of grid points used for the trapezoidal
#'   integration (default 200).
#' @return a one-row data frame with columns `est` and, if `se.fit =
#'   TRUE`, `se`, `lower`, `upper`.
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   time = rexp(200, 0.2), status = rbinom(200, 1, 0.7), x1 = rbinom(200, 1, 0.5)
#' )
#' fit <- rpsurv(survival::Surv(time, status) ~ x1, data = dat, tve = "x1")
#' rmst_diff(fit, newdata = data.frame(x1 = 1), newdata0 = data.frame(x1 = 0),
#'           tau = 4, se.fit = TRUE)
#' }
#' @export
rmst_diff <- function(object, newdata, newdata0, tau, se.fit = FALSE, n_grid = 200) {
  stopifnot(inherits(object, "rpsurv"), nrow(newdata) == 1, nrow(newdata0) == 1)
  grid <- seq(0, tau, length.out = n_grid)
  grid[1] <- max(grid[1], .Machine$double.eps) # avoid log(0) in the design
  log_grid <- log(grid)

  cov1 <- newdata[1, colnames(object$cov_data), drop = FALSE]
  cov0 <- newdata0[1, colnames(object$cov_data), drop = FALSE]
  d1 <- rp_predict_design(object, cov1, log_grid)
  d0 <- rp_predict_design(object, cov0, log_grid)
  sdiff <- d1$S - d0$S

  w <- diff(grid)
  trapw <- c(0, w / 2) + c(w / 2, 0) # standard trapezoidal node weights
  est <- sum(trapw * sdiff)

  out <- data.frame(est = est)
  if (se.fit) {
    grad <- rp_grad_surv(d1, object$scale) - rp_grad_surv(d0, object$scale)
    grad_w <- as.numeric(trapw %*% grad) # linear combination -> single gradient vector
    se <- sqrt(as.numeric(grad_w %*% object$vcov %*% grad_w))
    z <- stats::qnorm(0.975)
    out$se <- se
    out$lower <- est - z * se
    out$upper <- est + z * se
  }
  out
}
