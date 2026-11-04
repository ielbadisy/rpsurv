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
  design <- rp_design(log_times, cov_rep, object$knots, object$tve, object$tve_knots,
                      object$nle, object$nle_knots)
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

# Quadratic-form SE, one value per row of `grad` (n_times x p matrix). A
# negative quadratic form (possible only with a non-positive-definite vcov)
# gives NaN without R's generic "NaNs produced"; rp_check_output() then
# explains it.
#' @keywords internal
rp_delta_se <- function(grad, V) {
  q <- rowSums((grad %*% V) * grad)
  q[!is.na(q) & q < 0] <- NaN
  sqrt(q)
}

# Cumulative hazard -log S(eta), computed without underflow of S.
#' @keywords internal
rp_cumhaz <- function(eta, scale) {
  switch(scale,
    hazard = exp(eta),
    odds = ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))),
    normal = -stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
  )
}

#' @keywords internal
rp_check_times <- function(times, arg = "times") {
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times <= 0)) {
    stop("`", arg, "` must be positive, finite numbers (the model is defined on log time)", call. = FALSE)
  }
  invisible(times)
}

# Warn once, with the reason, when a prediction or effect curve contains
# non-finite values. `deta` holds d(eta)/d(log t) for hazard-based outputs:
# a non-positive value means the fitted hazard is not positive there.
#' @keywords internal
rp_check_output <- function(est, se = NULL, V = NULL, deta = NULL) {
  n_haz <- if (is.null(deta)) 0L else sum(!is.na(deta) & deta <= 0)
  if (n_haz) {
    warning(sprintf(
      "the fitted hazard is not positive at %d of %d evaluation point(s): the fitted cumulative hazard decreases there, typically when extrapolating beyond the observed follow-up. Hazard-based estimates there are negative, NaN or infinite.",
      n_haz, length(deta)
    ), call. = FALSE)
  }
  bad_est <- sum(is.nan(est) | is.infinite(est))
  if (bad_est && !n_haz) {
    warning(sprintf("%d estimate(s) are NaN or infinite (numerical overflow or underflow at extreme times or covariate values)",
                    bad_est), call. = FALSE)
  }
  if (!is.null(se)) {
    if (!is.null(V) && any(!is.finite(V))) {
      warning("the variance-covariance matrix is unavailable (see the warning from rpsurv()); standard errors and confidence limits are NA",
              call. = FALSE)
    } else {
      bad_se <- sum(is.nan(se) | is.infinite(se))
      if (bad_se && !n_haz) {
        warning(sprintf(
          "%d standard error(s) are NaN or infinite%s",
          bad_se,
          if (!is.null(V) && min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < 0) {
            " because the variance-covariance matrix is not positive definite (see the warning from rpsurv())"
          } else ""
        ), call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

# Covariate design columns (object$cov_data layout) for new data. `newdata`
# is normally given on the original variables of the formula (factors or
# character levels, untransformed numeric variables) and is coded through
# the fitted terms, factor levels and contrasts, as in predict.lm(). For
# backward compatibility, data that already holds the design columns (for
# example `celltypeadeno`) is used as is.
#' @keywords internal
rp_newdata_cov <- function(object, newdata, arg = "newdata") {
  newdata <- as.data.frame(newdata)
  cn <- colnames(object$cov_data)
  if (!length(cn)) return(newdata[, character(0), drop = FALSE])
  tt <- object$cov_terms
  vars <- if (is.null(tt)) cn else all.vars(tt)
  if (!is.null(tt) && all(vars %in% names(newdata))) {
    for (v in intersect(names(object$xlevels), names(newdata))) {
      x <- newdata[[v]]
      bad <- setdiff(unique(as.character(x[!is.na(x)])), object$xlevels[[v]])
      if (length(bad)) {
        stop(sprintf("`%s` has level(s) of %s not seen when fitting: %s (fitted levels: %s)",
                     arg, v, paste(bad, collapse = ", "),
                     paste(object$xlevels[[v]], collapse = ", ")), call. = FALSE)
      }
    }
    mm <- tryCatch({
      mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = object$xlevels)
      if (!is.null(cl <- attr(tt, "dataClasses"))) stats::.checkMFClasses(cl, mf)
      stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    }, error = function(e) {
      stop("`", arg, "` could not be coded like the fitted data: ", conditionMessage(e), call. = FALSE)
    })
    if (!all(cn %in% colnames(mm)) || nrow(mm) != nrow(newdata)) {
      stop("`", arg, "` could not be coded like the fitted data (design columns ",
           paste(setdiff(cn, colnames(mm)), collapse = ", "), " missing)", call. = FALSE)
    }
    return(as.data.frame(mm[, cn, drop = FALSE], optional = TRUE))
  }
  if (all(cn %in% names(newdata))) return(newdata[, cn, drop = FALSE])
  stop("`", arg, "` is missing covariate(s): ",
       paste(setdiff(vars, names(newdata)), collapse = ", "), call. = FALSE)
}

#' Predict from a fitted Royston-Parmar model
#'
#' @param object a fitted `"rpsurv"` object.
#' @param newdata a data frame of covariate values, one row per subject,
#'   on the original variables of the model formula: factors (or character
#'   values) by their levels, numeric variables untransformed, as for
#'   [survival::coxph()]. Factor levels not seen when fitting give an
#'   error. A data frame holding the design columns instead (for example
#'   `celltypeadeno` for the treatment-contrast indicator of a factor
#'   `celltype`) is also accepted.
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
  rp_check_times(times)
  log_times <- log(times)
  inv_link <- rp_inv_link(object$scale)
  V <- object$vcov
  all_deta <- list()
  all_se <- list()

  if (type %in% c("hr", "sdiff")) {
    if (is.null(newdata) || is.null(newdata0)) {
      stop("`type = \"", type, "\"` requires both `newdata` (contrast profile) and ",
           "`newdata0` (reference profile).", call. = FALSE)
    }
    newdata <- rp_newdata_cov(object, newdata, "newdata")
    newdata0 <- rp_newdata_cov(object, newdata0, "newdata0")
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
        all_deta[[i]] <- c(d1$deta, d0$deta)
        out <- data.frame(id = i, time = times, est = est)
        if (se.fit) {
          grad <- rp_grad_loghaz(d1, object$scale) - rp_grad_loghaz(d0, object$scale)
          se_log <- rp_delta_se(grad, V)
          all_se[[i]] <- se_log
          out$lower <- est * exp(-z * se_log)
          out$upper <- est * exp(z * se_log)
        }
      } else {
        est <- d1$S - d0$S
        out <- data.frame(id = i, time = times, est = est)
        if (se.fit) {
          grad <- rp_grad_surv(d1, object$scale) - rp_grad_surv(d0, object$scale)
          se_diff <- rp_delta_se(grad, V)
          all_se[[i]] <- se_diff
          out$lower <- est - z * se_diff
          out$upper <- est + z * se_diff
        }
      }
      results[[i]] <- out
    }
    res <- do.call(rbind, results)
    rp_check_output(res$est, if (se.fit) unlist(all_se), V,
                    if (type == "hr") unlist(all_deta))
    return(res)
  }

  newdata <- if (is.null(newdata)) object$cov_data else rp_newdata_cov(object, newdata)

  results <- vector("list", nrow(newdata))
  for (i in seq_len(nrow(newdata))) {
    cov_i <- newdata[i, colnames(object$cov_data), drop = FALSE]
    d <- rp_predict_design(object, cov_i, log_times)
    X <- d$X; eta <- d$eta

    est <- switch(type,
      link = eta,
      survival = d$S,
      cumhaz = rp_cumhaz(eta, object$scale),
      hazard = rp_hazard_from_design(d, times, object$scale)
    )
    if (type == "hazard") all_deta[[i]] <- d$deta

    out <- data.frame(id = i, time = times, est = est)

    if (se.fit && type %in% c("survival", "cumhaz", "link")) {
      se_eta <- rp_delta_se(X, V)
      all_se[[i]] <- se_eta
      lo_eta <- eta - z * se_eta
      hi_eta <- eta + z * se_eta
      if (type == "link") {
        out$lower <- lo_eta
        out$upper <- hi_eta
      } else if (type == "survival") {
        out$lower <- inv_link(hi_eta)
        out$upper <- inv_link(lo_eta)
      } else {
        # the cumulative hazard increases with eta
        out$lower <- rp_cumhaz(lo_eta, object$scale)
        out$upper <- rp_cumhaz(hi_eta, object$scale)
      }
    } else if (se.fit && type == "hazard") {
      grad <- rp_grad_loghaz(d, object$scale)
      se_log <- rp_delta_se(grad, V)
      all_se[[i]] <- se_log
      out$lower <- est * exp(-z * se_log)
      out$upper <- est * exp(z * se_log)
    }
    results[[i]] <- out
  }

  res <- do.call(rbind, results)
  rp_check_output(res$est, if (se.fit) unlist(all_se), V,
                  if (type == "hazard") unlist(all_deta))
  res
}

#' Restricted mean survival time difference between two covariate profiles
#'
#' Integrates the survival difference `sdiff = S(t | newdata) -
#' S(t | newdata0)` from 0 to `tau` by the trapezoidal rule, with an
#' optional delta-method standard error and 95% confidence interval.
#'
#' @param object a fitted `"rpsurv"` object.
#' @param newdata a one-row data frame, the contrast covariate profile
#'   (on the original variables, as in [predict.rpsurv()]).
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
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("`tau` must be a single positive, finite number", call. = FALSE)
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || n_grid < 2) {
    stop("`n_grid` must be a single number >= 2", call. = FALSE)
  }
  grid <- seq(0, tau, length.out = n_grid)
  grid[1] <- max(grid[1], .Machine$double.eps) # avoid log(0) in the design
  log_grid <- log(grid)

  cov1 <- rp_newdata_cov(object, newdata, "newdata")[1, , drop = FALSE]
  cov0 <- rp_newdata_cov(object, newdata0, "newdata0")[1, , drop = FALSE]
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
    se <- rp_delta_se(matrix(grad_w, nrow = 1L), object$vcov)
    z <- stats::qnorm(0.975)
    out$se <- se
    out$lower <- est - z * se
    out$upper <- est + z * se
  }
  rp_check_output(est, if (se.fit) out$se, object$vcov)
  out
}
