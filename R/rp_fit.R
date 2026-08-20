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

rp_start_values <- function(time, status, log_time, p, p_base, scale) {
  beta0 <- rep(0, p)
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
#' @param formula a `Surv(time, status) ~ covariates` formula.
#' @param data a data frame.
#' @param df degrees of freedom for the baseline spline in log time
#'   (number of interior knots is `df - 1`). Default 4.
#' @param knots optional numeric vector of ALL baseline knots (boundary
#'   knots first/last), overriding `df`.
#' @param tvc character vector of covariate names allowed a time-varying
#'   (non-proportional) effect, each modelled with its own spline in log time.
#' @param tvc.df degrees of freedom for each time-varying effect spline. Default 3.
#' @param scale one of `"hazard"` (proportional hazards, the Royston-Parmar
#'   default), `"odds"` (proportional odds) or `"normal"` (probit).
#' @param control list of control parameters passed to [stats::optim()].
#' @return an object of class `"rpsurv"`.
#' @export
rpsurv <- function(formula, data, df = 4, knots = NULL, tvc = NULL, tvc.df = 3,
                    scale = c("hazard", "odds", "normal"), control = list()) {
  scale <- match.arg(scale)
  scale_code <- rp_scale_code(scale)

  mf <- stats::model.frame(formula, data)
  resp <- stats::model.response(mf)
  if (!inherits(resp, "Surv")) stop("LHS of `formula` must be a Surv() object", call. = FALSE)
  if (ncol(resp) != 2L) stop("only right-censored data (Surv(time, status)) is currently supported", call. = FALSE)
  time <- resp[, 1L]
  status <- resp[, 2L]
  if (any(time <= 0)) stop("all survival times must be > 0", call. = FALSE)
  log_time <- log(time)

  cov_terms <- attr(stats::terms(formula), "term.labels")
  if (length(cov_terms)) {
    cov_data <- as.data.frame(stats::model.matrix(stats::reformulate(cov_terms, intercept = FALSE), mf))
  } else {
    cov_data <- data.frame(row.names = seq_along(time))[FALSE, , drop = FALSE]
    cov_data <- cov_data[seq_along(time), , drop = FALSE]
  }

  if (!is.null(tvc)) {
    missing_tvc <- setdiff(tvc, names(cov_data))
    if (length(missing_tvc)) {
      stop("`tvc` variable(s) not found among covariates: ", paste(missing_tvc, collapse = ", "), call. = FALSE)
    }
  }

  if (is.null(knots)) {
    knots <- default_knots(log_time[status == 1], df)
  } else {
    knots <- sort(knots)
  }

  tvc_knots <- NULL
  if (!is.null(tvc) && length(tvc)) {
    tvc_knots <- stats::setNames(
      lapply(tvc, function(v) default_knots(log_time[status == 1], tvc.df)),
      tvc
    )
  }

  design <- rp_design(log_time, cov_data, knots, tvc, tvc_knots)
  X <- design$X
  dX <- design$dX
  p <- ncol(X)
  p_base <- design$p_base

  beta0 <- rp_start_values(time, status, log_time, p, p_base, scale)

  negloglik <- function(beta) rp_negloglik_grad_cpp(beta, X, dX, log_time, as.numeric(status), scale_code)$value
  gradient  <- function(beta) rp_negloglik_grad_cpp(beta, X, dX, log_time, as.numeric(status), scale_code)$gradient

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
      tvc = tvc,
      tvc_knots = tvc_knots,
      formula = formula,
      call = match.call(),
      terms = stats::terms(formula),
      p_base = p_base,
      time = time,
      status = status,
      log_time = log_time,
      cov_data = cov_data,
      X = X,
      dX = dX
    ),
    class = "rpsurv"
  )
}
