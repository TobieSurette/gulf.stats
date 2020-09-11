#' Skewness
#'
#' @description Calculate sample skewness.
#'
#' @param x Numeric vector.
#' @param na.rm Logical value specifying whether to remove NA values from the calculation.
#'
#' @examples
#' x <- exp(rnorm(100))
#' skewness(x)

#' @export skewness
skewness <- function(x, na.rm = FALSE, ...){
   # Remove NA values:
   if (na.rm) x <- x[!base::is.na(x)]

   # Number of observations:
   n <- length(x)

   # Degenerate inputs:
   if (n < 3) return(NA)
   s <- stats::sd(x)
   if (s == 0) return(0)

   # Calculate sample skewness:
   v <- n /((n-1)*(n-2)) * sum(((x-mean(x))/s)^3)

   return(v)
}
