#' Covariance Functions
#'
#' @description Evaluate covariance functions.
#'
#' @param x \code{variogram} object.
#' @param h Lag distances.
#'
#'

#' @describeIn covariance Generic \code{variogram} method.
#' @export
covariance <- function(x, ...) UseMethod("covariance")

#' @describeIn covariance Evaluate covariance associated with variogram at specified distances.
#' @export
covariance.variogram <- function(x, h) v <- x$sill - variogram(x, h)

