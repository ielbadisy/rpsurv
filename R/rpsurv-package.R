#' rpsurv: Fast Royston-Parmar Flexible Parametric Survival Models
#'
#' @description
#' Fits Royston-Parmar flexible parametric survival models with the
#' likelihood and gradient evaluated in parallel C++ ('RcppParallel').
#'
#' @useDynLib rpsurv, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @keywords internal
"_PACKAGE"
