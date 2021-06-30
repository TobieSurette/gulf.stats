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
integrate.spline <- function(x, n = 1){
   if (n == 0) return(x)
   p <- attr(x, "polynomial")
   knots <- c(-Inf, attr(x, "knots"), Inf)
   for (i in 1:n){
      p[[1]] <- integrate(p[[1]])
      k <- p[[1]](knots[2]) - p[[1]](knots[1])
      for (j in 2:length(p)){
         p[[j]] <- integrate(p[[j]])
         p[[j]] <- p[[j]] + (-p[[j]](knots[j])) + k
         k <- p[[j]](knots[j+1])
      }
   }

   # Build spline object:
   f <- spline.default(polynomial = p, knots = attr(x, "knots"))

   return(f)
}
