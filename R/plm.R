#' @title Piecewise-Linear Model Class
#'
#' @description Set of functions which define a class of smooth piecewise linear models.
#'
#' @param x Numerical vector of values at which the function is to be evaluated.
#' @param scale Numeric value specifying the total width of the range of the basis functions.
#' @param skewness Degree of asymmetry in the basis functions. Negative values yield left-skewed basis functions,
#'                 right values yield right-skewed basis functions, while a zero value indicates no skewness.
#' @param degree Integer specifying the degree of the polynomial basis function.
#'
#' @examples
#'
#'

#' @export plm.basis
plm.basis <- function(x, scale = 1, skewness = 0, degree = 3, ...){
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

#' @export
plm <- function(x, knots = 0, intercept = 0, slope = c(0, 1), scale = 1, skewness = 0, ...){
   # Number of knots:
   k <- length(knots)

   # Parse input arguments:
   if (length(slope) != (k+1)) stop("'slope' must have one element more than 'knots'.")
   if (length(scale) == 1)    scale    <- rep(scale, k)
   if (length(skewness) == 1) skewness <- rep(skewness, k)
   if (length(scale) != k) stop("'scale' elemnts do no match the number of 'knots'.")
   if (length(skewness) != k) stop("'skewness' elemnts do no match the number of 'knots'.")

   # Define 'plm' function:
   fun <- function(x){
      y <- intercept + slope[1] * x
      if (length(k > 0)){
         # Generate basis functions:
         for (i in 1:k){
            b <- plm.basis(scale = scale[i], slope = slope[i], skewness = skewness[i])
            b <- integrate.spline(b, n = 2)
            y <- y + (slope[i+1]-slope[i]) * b(x - knots[i])
         }
      }
      return(y)
   }

   # Return result:
   if (missing(x)) return(fun) else return(fun(x))
}
