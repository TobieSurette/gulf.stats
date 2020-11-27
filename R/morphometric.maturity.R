#' Determine Morphometric Maturity
#'
#' @description Functions to determine sexual maturity from morphometric measurements.
#'
#' @param x Data object.
#' @param probability Logical value specifying whether the morphoemtric data are rounded.
#'
#' @examples
#' x <- read.scsbio(2010)
#' v <- morphometric.maturity(x)
#'
#' @seealso \code{\link{morphometry}}
#' @seealso \code{\link[gulf.data]{maturity}}

#' @describeIn morphometric.maturity Generic morphometric maturity method.
#' @export morphometric.maturity
morphometric.maturity <- function(x, ...) UseMethod("morphometric.maturity")

#' @describeIn morphometric.maturity Determine morphometric maturity from snow crab biological data.
#' @export morphometric.maturity
morphometric.maturity.scsbio <- function(x, probability = FALSE, ...){
   # Survey years:
   years <- sort(unique(year(x)))

   v <- rep(NA, nrow(x))
   for (i in 1:2){
      for (j in 1:length(years)){
         # Read data:
         ix <- which((year(x) == years[j]) & (x$sex == i) & !is.na(b$carapace.width))

         if (length(ix)){
            # Prepare variables:
            if (i == 1) y <- x$chela.height[ix]
            if (i == 2) y <- x$abdomen.width[ix]

            # Unconditional maturity identification:
            z <- rep(NA, nrow(x))
            z[which(x < 30)] <- 0 # Crab smaller than 30mm CW are considered immature.

            # Fit morphometric model:
            theta <- fit.morphometry.scsbio(x$carapace.width[ix], y, z, sex = i, discrete = years[j] < 1998)

            # Extract maturity proportions:
            v[ix] <- morphometry.scsbio(x, y, theta = theta, discrete = years[j] < 1998)$p_mature_posterior
         }
      }
   }

   # Convert to binary maturity:
   if (!probability) v <- (v >= 0.5)

   return(v)
}
