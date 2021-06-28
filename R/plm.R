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

plm <- function(x, scale = 1, skewness = 0, degree = 3, ...){

  return(NULL)
}
