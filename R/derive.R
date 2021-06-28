#' @title Derivation
#' 
#' @description Methods to perform analytic or numeric derivation.

#' @export
derive <- function(x) UseMethod("derive")

#' @export 
#' @describeIn polynomial Polynomial derivation method.
derive.polynomial <- function(x){  
  beta <- coef(x)
  if (length(beta) == 1) return(polynomial(0))
  beta <- beta[-1] * (1:(length(beta)-1) )
  p <- polynomial(beta)
  return(p)
}

#' @export 
#' @describeIn polynomial Spline derivation method.
derive.spline <- function(x){  
  p <- attr(x, "polynomial")
  for (i in 1:length(p)) p[[i]] <- derive(p[[i]])
  f <- spline.default(polynomial = p, knots = attr(x, "knots"))
  
  return(f)
}
