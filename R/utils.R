#' Null coalescing operator
#'
#' Returns the left-hand side if not NULL, otherwise returns the right-hand side.
#'
#' @name or_or
#' @aliases %||%
#' @param a Left side value
#' @param b Right side value (default)
#' @return a if not NULL, otherwise b
#' @keywords internal
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a
