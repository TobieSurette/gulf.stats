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

#' @export morphometric.maturity
morphometric.maturity <- function(x, ...) UseMethod("morphometric.maturity")

#' @describeIn morphometric.maturity Determine morphometric maturity from snow crab biological data.
#' @rawNamespace S3method(morphometric.maturity,scsbio)
morphometric.maturity.scsbio <- function(x, probability = FALSE, ...){
   # Survey years:
   years <- sort(unique(year(x)))

   v <- rep(NA, nrow(x))
   for (i in 1:2){
      for (j in 1:length(years)){
         # Read data:
         ix <- which((year(x) == years[j]) & (x$sex == i) & !is.na(x$carapace.width))

         if (length(ix) > 0){
            # Prepare variables:
            if (i == 1) y <- x$chela.height[ix]
            if (i == 2) y <- x$abdomen.width[ix]

            # Round off to 0.2mm to increase optimization speed:
            r <- aggregate(list(n = y), list(x = round(2*x$carapace.width[ix])/2, y = round(2*y)/2), length)

            # Unconditional maturity identification:
            z <- rep(NA, length(ix))
            z[which(x$carapace.width[ix] < 30)] <- 0 # Crab smaller than 30mm CW are considered immature.

            # Fit morphometric model:
            theta <- fit.morphometry.scsbio(x$carapace.width[ix], y, z, n, sex = i, discrete = years[j] < 1998)

            # Extract maturity proportions:
            v[ix] <- morphometry(x$carapace.width[ix], y, z = z, theta = theta, discrete = years[j] < 1998)$p_mature_posterior
         }
      }
   }

   # Convert to binary maturity:
   if (!probability) v <- (v >= 0.5)

   return(v)
}
