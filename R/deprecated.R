#' Deprecated functions in rpsurv
#'
#' These functions were renamed to glued lowercase names and will be
#' removed in a future release. Use the new names instead.
#'
#' @name rpsurv-deprecated
#' @param ... Arguments passed to the new function.
#' @return The value of the new function.
#' @keywords internal
NULL

#' @rdname rpsurv-deprecated
#' @export
coxsnell_plot <- function(...) {
  .Deprecated("coxsnellplot")
  coxsnellplot(...)
}

#' @rdname rpsurv-deprecated
#' @export
km_compare_plot <- function(...) {
  .Deprecated("kmcompareplot")
  kmcompareplot(...)
}

#' @rdname rpsurv-deprecated
#' @export
rmst_diff <- function(...) {
  .Deprecated("rmstdiff")
  rmstdiff(...)
}
