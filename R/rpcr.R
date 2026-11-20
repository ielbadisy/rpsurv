#' Flexible parametric Fine-Gray model for competing risks
#'
#' Fits a Royston-Parmar model for the subdistribution hazard of one cause (Fine and Gray, 1999), so
#' that the cumulative incidence function (CIF) is smooth and modelled directly:
#' \deqn{\log\{-\log(1 - F_k(t \mid x))\} = s(\log t; \gamma) + x^\top \beta,}
#' with `exp(beta)` the subdistribution hazard ratios (sHR). As in Lambert, Wilkes and Crowther
#' (2017), subjects with a competing event stay in the risk set after it, weighted by the inverse
#' probability of remaining uncensored, and the weighted (counting-process) Royston-Parmar likelihood
#' is maximised. Unlike [survival::finegray()], whose expansion is designed for the Cox partial
#' likelihood, the data are split at every censoring time, so that the weights are exact over the
#' whole follow-up. The variance is the robust sandwich clustered on the subject, which agrees closely
#' with that of `cmprsk::crr()`.
#'
#' @param formula `Surv(time, status) ~ covariates`, where `status` is a factor whose first level is
#'   censoring and whose other levels are the event types, as for [survival::finegray()] (for a
#'   numeric code 0/1/2, use `factor(status, 0:2)`).
#' @param data a data frame, one row per subject.
#' @param cause the event type (a level of `status`) whose cumulative incidence is modelled.
#' @param df,knots,tve,tve.df,nle,nle.df,scale,control as in [rpsurv()], applied to the
#'   subdistribution: `scale = "hazard"` is the Fine-Gray (proportional subdistribution hazards)
#'   model, `"odds"` gives proportional odds of the CIF, and `tve` lets a subdistribution hazard ratio
#'   change over time. The default knots are placed at centiles of the log event times of `cause`.
#' @details The fit is an `"rpsurv"` object for the subdistribution "survival" `1 - F_k(t)`, so
#'   [coef()], [vcov()] (robust), [confint()], [summary()], [nlecurve()] and [tvecurve()] apply as
#'   usual. Its log-likelihood is a weighted pseudo-likelihood: [AIC()] and [BIC()] compare models of
#'   the same cause on the same data only. The naive (model-based) variance is kept in `naive.vcov`.
#' @return an object of class `c("rpcr", "rpsurv")`, with in addition `cause` and `naive.vcov`.
#' @references Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a
#'   competing risk. JASA 94:496-509.
#'
#'   Lambert PC, Wilkes SR, Crowther MJ (2017). Flexible parametric modelling of the cause-specific
#'   cumulative incidence function. Statistics in Medicine 36:1429-1446.
#' @seealso [predict.rpcr()], [plot.rpcr()].
#' @examples
#' \donttest{
#' d <- survival::mgus2
#' d$etime <- with(d, ifelse(pstat == 1, ptime, futime))
#' d$event <- factor(with(d, ifelse(pstat == 1, 1, 2 * death)), 0:2, c("censor", "pcm", "death"))
#' fit <- rpcr(survival::Surv(etime, event) ~ sex + age, data = d, cause = "pcm")
#' summary(fit)
#' predict(fit, data.frame(sex = c("F", "M"), age = 70), times = c(120, 240), se.fit = TRUE)
#' plot(fit, newdata = data.frame(sex = c("F", "M"), age = 70))
#' }
#' @export
rpcr <- function(formula, data, cause = 1, df = 4, knots = NULL, tve = NULL, tve.df = 3,
                 nle = NULL, nle.df = 3, scale = c("hazard", "odds", "normal"), control = list()) {
  scale <- match.arg(scale)
  labs <- attr(stats::terms(formula, data = data), "term.labels")
  mf <- stats::model.frame(formula, data)
  Y <- stats::model.response(mf)
  if (!inherits(Y, "Surv") || attr(Y, "type") != "mright") {
    stop("the response must be Surv(time, status) with `status` a factor (first level censoring)",
         call. = FALSE)
  }
  k <- match(as.character(cause), attr(Y, "states"))
  if (is.na(k)) stop("`cause` must be one of: ", paste(attr(Y, "states"), collapse = ", "), call. = FALSE)
  om <- as.integer(attr(mf, "na.action"))
  sub <- if (length(om)) data[-om, , drop = FALSE] else data
  e <- rp_fgexpand(Y[, 1L], Y[, 2L], k)
  fg <- cbind(sub[e$id, , drop = FALSE], .fgstart = e$start, .fgstop = e$stop, .fgstatus = e$status)
  f <- stats::reformulate(if (length(labs)) labs else "1",
                          quote(survival::Surv(.fgstart, .fgstop, .fgstatus)), env = environment(formula))
  fit <- rpsurv(f, fg, df = df, knots = knots, tve = tve, tve.df = tve.df, nle = nle, nle.df = nle.df,
                scale = scale, control = control, weights = e$w)
  # Sandwich clustered on the subject: a subject with a competing event contributes several weighted
  # rows, and the IPCW weights make the model-based variance invalid.
  U <- rowsum(fit$weights * rp_scores(fit), e$id)
  V <- fit$vcov %*% crossprod(U) %*% fit$vcov
  fit[c("naive.vcov", "vcov", "cause", "n", "call", "formula")] <-
    list(fit$vcov, (V + t(V)) / 2, as.character(cause), nrow(sub), match.call(), formula)
  fit$na.action <- attr(mf, "na.action")
  fit$counting <- FALSE
  fit$cov_data <- rp_newdata_cov(fit, sub)
  fit$cov_vars <- if (!is.null(fit$cov_terms)) stats::get_all_vars(fit$cov_terms, sub)
  class(fit) <- c("rpcr", class(fit))
  fit
}

# Fine-Gray data expansion. A subject with a competing event at T stays at risk after T, with weight
# G(t) / G(T) (G the Kaplan-Meier of censoring), up to the last follow-up time. The rows are split at
# every censoring time, where G steps, so the weight is exact over each row: survival::finegray()
# splits only at event times of the cause, which is enough for the Cox partial likelihood but not
# for the cumulative hazard integral of a parametric likelihood. `ev` is 0 for censoring, else the
# event type; returns the subject index, (start, stop], status and weight of each row.
#' @keywords internal
rp_fgexpand <- function(time, ev, cause) {
  sf <- survival::survfit(survival::Surv(time, ev == 0) ~ 1)
  ct <- sf$time[sf$n.event > 0]
  G <- stats::stepfun(ct, c(1, sf$surv[sf$n.event > 0]))
  tmax <- max(time)
  cr <- which(ev != 0 & ev != cause & time < tmax)
  ext <- lapply(cr, function(i) {
    s <- c(time[i], ct[ct > time[i] & ct < tmax])
    cbind(i, s, c(s[-1L], tmax), G(s) / G(time[i]))
  })
  ext <- do.call(rbind, c(list(matrix(0, 0, 4)), ext))
  ext <- ext[ext[, 4L] > 0, , drop = FALSE]
  n <- length(time)
  o <- order(c(seq_len(n), ext[, 1L]))
  list(id = c(seq_len(n), ext[, 1L])[o], start = c(rep(0, n), ext[, 2L])[o],
       stop = c(time, ext[, 3L])[o], status = c(as.numeric(ev == cause), rep(0, nrow(ext)))[o],
       w = c(rep(1, n), ext[, 4L])[o])
}

# Per-row score (gradient of the unweighted log-likelihood term wrt beta), the row analogue of
# rp_negloglik_grad_cpp().
#' @keywords internal
rp_scores <- function(object) {
  b <- object$coefficients
  sc <- function(eta) switch(object$scale,
    hazard = list(dl = -exp(eta), gp = 1),
    odds = list(dl = -stats::plogis(eta), gp = stats::plogis(-eta)),
    normal = {
      r <- exp(stats::dnorm(eta, log = TRUE) - stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE))
      list(dl = -r, gp = r - eta)
    })
  X <- object$X
  s <- sc(drop(X %*% b))
  e <- sc(drop(object$Xentry %*% b))$dl * object$has_entry
  object$status * (object$dX / drop(object$dX %*% b) + s$gp * X) + s$dl * X - e * object$Xentry
}

#' Predict from a flexible parametric Fine-Gray model
#'
#' @param object a fitted `"rpcr"` object.
#' @param newdata,newdata0,times as in [predict.rpsurv()]. `newdata` defaults to the subjects used in
#'   the fit (one curve per subject).
#' @param type one of `"cif"` (cumulative incidence of the cause), `"cifdiff"` (CIF difference
#'   `F(t | newdata) - F(t | newdata0)`), `"shr"` (subdistribution hazard ratio of `newdata` vs
#'   `newdata0`, time dependent under `tve`) or `"link"` (the linear predictor).
#' @param se.fit logical; if `TRUE`, adds 95% delta-method limits `lower` and `upper` from the robust
#'   variance: on the linear-predictor scale for `"cif"` (so they stay in `[0, 1]`), the log scale for
#'   `"shr"` and the natural scale for `"cifdiff"`.
#' @param ... unused.
#' @return a long data frame with columns `id` (row of `newdata`), `time`, `est` and optionally
#'   `lower`, `upper`.
#' @export
predict.rpcr <- function(object, newdata = NULL, newdata0 = NULL, times = NULL,
                         type = c("cif", "cifdiff", "shr", "link"), se.fit = FALSE, ...) {
  type <- match.arg(type)
  p <- predict.rpsurv(object, newdata, newdata0, times,
                      switch(type, cif = "survival", cifdiff = "sdiff", shr = "hr", link = "link"), se.fit)
  if (type %in% c("cif", "cifdiff")) {
    a <- as.numeric(type == "cif")
    p$est <- a - p$est
    if (se.fit) p[c("lower", "upper")] <- list(a - p$upper, a - p$lower)
  }
  p
}

#' Plot cumulative incidence curves from a flexible parametric Fine-Gray model
#'
#' @param x a fitted `"rpcr"` object.
#' @param newdata covariate profile(s), one row per curve; defaults to one profile with numeric
#'   variables at their mean and factors at their most common value, as in [plot.rpsurv()].
#' @param type `"cif"` (cumulative incidence) or `"shr"` (subdistribution hazard ratio of each row of
#'   `newdata` vs `newdata0`, useful under `tve`).
#' @param newdata0 reference profile(s) for `type = "shr"`.
#' @param times,ci,col,legend,labels,... as in [plot.rpsurv()]: `legend` is `TRUE`, `FALSE` or a
#'   position, and custom curve labels go in `labels = c("Remission", "Relapse")` (or `legend =` a
#'   vector of labels).
#' @return invisibly, the plotted predictions (see [predict.rpcr()]).
#' @export
plot.rpcr <- function(x, newdata = NULL, type = c("cif", "shr"), newdata0 = NULL, times = NULL,
                      ci = TRUE, col = NULL, legend = NULL, labels = NULL, ...) {
  type <- match.arg(type)
  newdata <- as.data.frame(if (is.null(newdata)) rp_default_profile(x) else newdata)
  n <- max(1L, nrow(newdata))
  col <- rep_len(if (is.null(col)) seq_len(n) else col, n)
  pred <- predict(x, newdata, newdata0, times, type = type, se.fit = ci)
  ylab <- if (type == "cif") paste0("Cumulative incidence (cause ", x$cause, ")") else
    "Subdistribution hazard ratio"
  rp_draw_curves(pred, newdata, ylab, type == "cif", col, legend, list(...), labels)
}
