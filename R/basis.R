#' @title Basis Functions
#'
#' @description Basis function definitions.
#'
#' @param x Predictor vector.
#' @param scale Scale parameter (e.g. width of the interval over which the basis function is defined.).
#' @param skewness Skewness parameter of basis function. Zero value implies no skewness, negative values are
#'                 left-skewed and positive values are right-skewed.
#' @param degree Degree of polynomial basis.
#

#' @export
basis <- function(x, ...) UseMethod("basis")

#' @export
basis.plm <- function(x, scale = 1, skewness = 0, degree = 3, ...){
   # Convert skewness to proportions:
   w <- 1 / (1+exp(-skewness))

   # Define knots:
   lower <- (w-1) * scale
   upper <- w * scale

   # Define basis:
   if (degree == 3){
      k <- (2 / (abs(lower) + abs(upper)))
      p <- c(polynomial(0),
             k * polynomial(c(1, 0, -3/lower^2, 2/lower^3)),
             k * polynomial(c(1, 0, -3/upper^2, 2/upper^3)),
             polynomial(0))
      knots <- c(lower, 0, upper)
      f <- spline.default(polynomial = p, knots = knots)
   }

   return(f)
}
