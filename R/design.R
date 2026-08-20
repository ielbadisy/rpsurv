#' Build the design matrices for a Royston-Parmar model
#'
#' Constructs `eta`'s design matrix `X` and the design matrix for
#' `d(eta)/d(log t)`, `dX`, which share the same column layout:
#' `[intercept, baseline spline(log t), covariates, tvc spline(log t):covariate]`.
#'
#' @keywords internal
rp_design <- function(log_time, cov_data, knots, tvc = NULL, tvc_knots = NULL) {
  p_base <- ncol(rcs_basis(log_time, knots))
  X_base <- cbind(1, rcs_basis(log_time, knots))
  dX_base <- cbind(0, rcs_basis(log_time, knots, derivative = TRUE))
  colnames(X_base) <- colnames(dX_base) <- c("(Intercept)", paste0("s(logt)", seq_len(p_base)))

  has_cov <- !is.null(cov_data) && ncol(cov_data) > 0L
  if (has_cov) {
    X_cov <- as.matrix(cov_data)
    dX_cov <- matrix(0, nrow(X_cov), ncol(X_cov))
    colnames(dX_cov) <- colnames(X_cov)
  } else {
    X_cov <- matrix(0, length(log_time), 0L)
    dX_cov <- X_cov
  }

  X_tvc <- NULL
  dX_tvc <- NULL
  if (!is.null(tvc) && length(tvc)) {
    X_tvc_list <- vector("list", length(tvc))
    dX_tvc_list <- vector("list", length(tvc))
    for (i in seq_along(tvc)) {
      var <- tvc[i]
      k <- tvc_knots[[var]]
      basis <- rcs_basis(log_time, k)
      dbasis <- rcs_basis(log_time, k, derivative = TRUE)
      x <- cov_data[[var]]
      X_tvc_list[[i]] <- basis * x
      dX_tvc_list[[i]] <- dbasis * x
      colnames(X_tvc_list[[i]]) <- colnames(dX_tvc_list[[i]]) <-
        paste0(var, ":s(logt)", seq_len(ncol(basis)))
    }
    X_tvc <- do.call(cbind, X_tvc_list)
    dX_tvc <- do.call(cbind, dX_tvc_list)
  }

  X <- cbind(X_base, X_cov, X_tvc)
  dX <- cbind(dX_base, dX_cov, dX_tvc)
  list(X = X, dX = dX, p_base = p_base + 1L)
}
