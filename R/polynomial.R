#' Polynomial evaluation
#' 
#' @description Horner's method for polynomial evaluation.
#' 
#' @param x Numeric vector at which the polynomial will be evaluated.
#' @param p Polynomial coefficients. The first coefficient corresponds to the intercept 
#'          or constant coefficent, while the last element corresponds to the coefficient
#'          of the highest polynomial order.
#'          
#' @examples 
#' x <- 1:10
#' polyval(x, 1) # Constant function.
#' polyval(x, c(1,0,2)) # Quadratic 1 + 2*x^2.
 
#' @export
polyval <- function(x, p){
   p <- rev(p) # Intercept is last parameter:
   y <- p[[1]] * rep(1, length(x))  
   if (length(p) > 1) for (i in 2:length(p)) y <- y * x + p[i]
   return(y)
} 
