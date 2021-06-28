#' @title Integration
#'
#' @description Methods to perform analytic or numeric integration.


#' @export
integrate <- function(x) UseMethod("integrate")

#' @export
integrate.default <- stats::integrate

#' @export
#' @describeIn polynomial Polynomial integration method.
integrate.polynomial <- function(x){
   beta <- c(0, coef(x) / (1:length(coef(x))))
   p <- polynomial(beta)
   return(p)
}

#' @export
#' @describeIn polynomial Spline integration method.
integrate.spline <- function(x){
   p <- attr(x, "polynomial")
   knots <- c(-Inf, attr(x, "knots"), Inf)
   k <- 0
   for (i in 1:length(p)){
      p[[i]] <- integrate(p[[i]])
      k <- k + p[[i]](knots[i+1]) - p[[i]](knots[i])
      p[[i]] <- k * p[[i]]
   }

   f <- spline.default(polynomial = p, knots = attr(x, "knots"))

   return(f)
}
