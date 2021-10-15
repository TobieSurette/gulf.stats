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

#' @export
plm <- function(x, theta, knots = 0, intercept = 0, slope = c(0, 1), scale = 1, skewness = 0, ...){
   # Parse parameter vector:
   if (!missing(theta)){
      theta[grep("^log_", names(theta))] <- exp(theta[grep("^log_", names(theta))])
      knots     <- theta[grep("knots*[0-9]", names(theta))]
      slope     <- theta[sort(grep("slopes*[0-9]", names(theta)))]
      scale     <- theta[sort(grep("scales*[0-9]+", names(theta)))]
      skewness  <- theta[sort(grep("skewness", names(theta)))]
      intercept <- theta[sort(grep("intercept", names(theta)))]
   }
   knots <- sort(knots)

   # Number of knots:
   k <- length(knots)

   # Parse input arguments:
   if (length(scale) == 0) scale <- 1
   if (length(skewness) == 0) skewness <- 0
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
            b <- basis.plm(scale = scale[i], slope = slope[i], skewness = skewness[i])
            b <- integrate.spline(b, n = 2)
            y <- y + (slope[i+1]-slope[i]) * b(x - knots[i])
         }
      }
      return(y)
   }

   # Return result:
   if (missing(x)) return(fun) else return(fun(x))
}

