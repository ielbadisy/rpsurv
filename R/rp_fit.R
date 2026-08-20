rp_scale_code <- function(scale) {
  switch(scale,
    hazard = 0L, ph = 0L,
    odds = 1L, po = 1L,
    normal = 2L, probit = 2L,
    stop("`scale` must be one of \"hazard\", \"odds\", \"normal\"", call. = FALSE)
  )
}

rp_link_fun <- function(scale) {
  switch(scale,
    hazard = , ph = function(S) log(-log(S)),
    odds = , po = function(S) log((1 - S) / S),
    normal = , probit = function(S) stats::qnorm(1 - S)
  )
}

rp_start_values <- function(time, status, p, scale, max_n = 20000L) {
  beta0 <- rep(0, p)
  if (length(time) > max_n) {
    idx <- as.integer(seq.int(1L, length(time), length.out = max_n))
    time <- time[idx]
    status <- status[idx]
  }
  sf <- survival::survfit(survival::Surv(time, status) ~ 1)
  keep <- sf$surv > 0 & sf$surv < 1 & sf$time > 0
  if (sum(keep) >= 2L) {
    link <- rp_link_fun(scale)
    y <- link(sf$surv[keep])
    x <- log(sf$time[keep])
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) >= 2L) {
      fit <- stats::lm(y[ok] ~ x[ok])
      beta0[1L] <- unname(stats::coef(fit)[1L])
      beta0[2L] <- unname(stats::coef(fit)[2L])
    }
  }
  if (!is.finite(beta0[2L]) || beta0[2L] <= 0) beta0[2L] <- 1
  if (!is.finite(beta0[1L])) beta0[1L] <- 0
  beta0
}

#' Fit a Royston-Parmar flexible parametric survival model
#'
#' @param formula a survival formula. Either `Surv(time, status) ~ covariates`
#'   for standard right-censored data, or `Surv(start, stop, status) ~ covariates`
#'   (counting-process form) for left-truncated data and/or genuine
#'   **time-varying covariates** (a covariate whose *value* changes over
#'   follow-up): split each subject's follow-up into intervals of constant
#'   covariate values and pass one row per interval. See Details.
#' @param data a data frame.
#' @param df degrees of freedom for the baseline spline in log time
#'   (number of interior knots is `df - 1`). Default 4.
#' @param knots optional numeric vector of ALL baseline knots (boundary
#'   knots first/last), overriding `df`.
#' @param tve character vector of covariate names allowed a **time-varying
#'   effect** (non-proportional hazards/odds/probit): the covariate's
#'   *coefficient* beta(t) is modelled as its own spline in log time, while
#'   its value is still fixed for a subject. This is distinct from a
#'   time-varying *covariate* (see `formula`); the two can be combined.
#' @param tve.df degrees of freedom for each time-varying effect spline. Default 3.
#' @param scale one of `"hazard"` (proportional hazards, the Royston-Parmar
#'   default), `"odds"` (proportional odds) or `"normal"` (probit).
#' @param control list of control parameters passed to [stats::optim()].
#'
#' @details
#' # Time-varying effect vs. time-varying covariate
#' These are two distinct extensions and `rpsurv()` supports both, separately
#' or combined:
#' * **Time-varying effect** (non-proportional hazards): a covariate's value
#'   is fixed for a subject but its *association with the outcome* changes
#'   over time, e.g. a treatment effect that fades. Request this with `tve`.
#' * **Time-varying covariate**: a covariate's *value* itself changes during
#'   follow-up, e.g. a lab measurement updated at clinic visits. Represent
#'   this in the data as counting-process (start, stop] intervals, one row
#'   per interval with the covariate value held constant within it, and fit
#'   with `Surv(start, stop, status) ~ ...`. `rpsurv()` handles the
#'   corresponding left-truncated likelihood automatically.
#' @return an object of class `"rpsurv"`.
#' @export
rpsurv <- function(formula, data, df = 4, knots = NULL, tve = NULL, tve.df = 3,
                    scale = c("hazard", "odds", "normal"), control = list()) {
  scale <- match.arg(scale)
  scale_code <- rp_scale_code(scale)

  mf <- stats::model.frame(formula, data)
  resp <- stats::model.response(mf)
  if (!inherits(resp, "Surv")) stop("LHS of `formula` must be a Surv() object", call. = FALSE)
  if (!ncol(resp) %in% c(2L, 3L)) {
    stop("`formula` must use Surv(time, status) or Surv(start, stop, status)", call. = FALSE)
  }
  counting <- ncol(resp) == 3L
  if (counting) {
    entry <- resp[, 1L]
    time <- resp[, 2L]
    status <- resp[, 3L]
  } else {
    entry <- rep(0, nrow(resp))
    time <- resp[, 1L]
    status <- resp[, 2L]
  }
  if (any(time <= 0)) stop("all exit times must be > 0", call. = FALSE)
  if (any(entry < 0) || any(entry >= time)) stop("entry times must satisfy 0 <= start < stop", call. = FALSE)
  log_time <- log(time)
  has_entry <- as.numeric(entry > 0)
  status_num <- as.numeric(status)
  log_entry_safe <- entry
  log_entry_safe[log_entry_safe <= 0] <- 1
  log_entry_safe <- log(log_entry_safe)

  cov_terms <- attr(stats::terms(formula), "term.labels")
  if (length(cov_terms)) {
    cov_data <- as.data.frame(stats::model.matrix(stats::reformulate(cov_terms, intercept = FALSE), mf))
  } else {
    cov_data <- data.frame(row.names = seq_along(time))[FALSE, , drop = FALSE]
    cov_data <- cov_data[seq_along(time), , drop = FALSE]
  }

  if (!is.null(tve)) {
    missing_tve <- setdiff(tve, names(cov_data))
    if (length(missing_tve)) {
      stop("`tve` variable(s) not found among covariates: ", paste(missing_tve, collapse = ", "), call. = FALSE)
    }
  }

  if (is.null(knots)) {
    knots <- default_knots(log_time[status == 1], df)
  } else {
    knots <- sort(knots)
  }

  tve_knots <- NULL
  if (!is.null(tve) && length(tve)) {
    tve_knots <- stats::setNames(
      lapply(tve, function(v) default_knots(log_time[status == 1], tve.df)),
      tve
    )
  }

  design <- rp_design(log_time, cov_data, knots, tve, tve_knots)
  X <- design$X
  dX <- design$dX
  p <- ncol(X)
  p_base <- design$p_base

  Xentry <- if (any(has_entry > 0)) {
    rp_design(log_entry_safe, cov_data, knots, tve, tve_knots)$X
  } else {
    matrix(0, nrow(X), p)
  }

  beta0 <- rp_start_values(time, status, p, scale)

  # optim() calls fn() then gr() at the same trial point on every BFGS step;
  # cache the (fused) C++ call so each point is only evaluated once
  cache <- new.env(parent = emptyenv())
  evaluate <- function(beta) {
    if (is.null(cache$beta) || !identical(beta, cache$beta)) {
      cache$result <- rp_negloglik_grad_cpp(beta, X, dX, Xentry, has_entry, log_time, status_num, scale_code)
      cache$beta <- beta
    }
    cache$result
  }
  negloglik <- function(beta) evaluate(beta)$value
  gradient  <- function(beta) evaluate(beta)$gradient

  opt_control <- utils::modifyList(list(maxit = 500, reltol = 1e-10), control)
  fit <- stats::optim(beta0, negloglik, gradient, method = "BFGS",
                       control = opt_control, hessian = TRUE)

  vcov <- tryCatch(solve(fit$hessian), error = function(e) {
    warning("Hessian not invertible; standard errors unavailable", call. = FALSE)
    matrix(NA_real_, p, p)
  })

  structure(
    list(
      coefficients = stats::setNames(fit$par, colnames(X)),
      vcov = `dimnames<-`(vcov, list(colnames(X), colnames(X))),
      loglik = -fit$value,
      n = length(time),
      nevent = sum(status),
      df = p,
      convergence = fit$convergence,
      scale = scale,
      knots = knots,
      tve = tve,
      tve_knots = tve_knots,
      counting = counting,
      formula = formula,
      call = match.call(),
      terms = stats::terms(formula),
      p_base = p_base,
      time = time,
      entry = entry,
      has_entry = has_entry,
      status = status,
      log_time = log_time,
      cov_data = cov_data,
      X = X,
      dX = dX,
      Xentry = Xentry
    ),
    class = "rpsurv"
  )
}
