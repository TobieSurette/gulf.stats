#' @title Calculate Error
#'
#' @description Calculate regression error for various statistical models.
#'
#' @param x Predictor vector.
#' @param theta Parameter vector.
#'

#' @export
sigma <- function(x, ...) UseMethod("sigma")

#' @export
sigma.plm <- function(x, theta, ...){
   if (missing(theta)) stop("'theta' must be specified.")

   # Define error function:
   fun <- function(x){
      # Calculate error function:
      v <- exp(theta[["log_sigma"]])
      if ("log_sigma_descent" %in% names(theta)){
         bound <- theta[["knots1"]] - exp(theta[["log_scale1"]])/2 - exp(theta[["log_scale_descent"]]) / 2
         b_descent <- integrate(basis.plm(scale = exp(theta[["log_scale_descent"]])))
         v <- v + exp(theta[["log_sigma_descent"]]) * (1-b_descent(x - bound))
      }
      if ("log_sigma_ascent" %in% names(theta)){
         bound <- theta[["knots2"]] + exp(theta[["log_scale2"]])/2 + exp(theta[["log_scale_ascent"]]) / 2
         b_ascent <- integrate(basis.plm(scale = exp(theta[["log_scale_ascent"]])))
         v <- v + exp(theta[["log_sigma_ascent"]]) * b_ascent(x - bound)
      }

      return(v)
   }

   # Return results:
   if (missing(x)) return(fun) else return(fun(x))
}
