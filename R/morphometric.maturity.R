#' Determine Morphometric Maturity
#'
#' @description Functions to determine sexual maturity from morphometric measurements.
#'
#' @param x Data object.
#' @param probability Logical value specifying whether the morphometric data are rounded.
#'
#' @examples
#' # Male morphometric maturity:
#' x <- read.scsbio(2011, sex = 1)
#' v <- morphometric.maturity.scsbio(x)
#' plot(c(0, 140), c(0, 40), type = "n", xlab = "Carapace width(m)", ylab = "Chela height(m)", xaxs = "i", yaxs = "i")
#' grid()
#' points(x$carapace.width[v], x$chela.height[v], pch = 21, bg = "red", cex = 0.2, lwd = 0.2)
#' points(x$carapace.width[!v], x$chela.height[!v], pch = 21, bg = "blue", cex = 0.2, lwd = 0.2)
#' box()
#'
#' @seealso \code{\link{morphometry}}
#' @seealso \code{\link[gulf.data]{maturity}}

#' @export morphometric.maturity
morphometric.maturity <- function(x, ...) UseMethod("morphometric.maturity")

#' @describeIn morphometric.maturity Determine morphometric maturity from snow crab biological data.
#' @export morphometric.maturity.scsbio
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
            bin <- 0.5
            y <- round(y / bin) * bin
            r <- aggregate(list(n = y), list(x = round(x$carapace.width[ix] / bin) * bin, y = y), length)
            tmp <- table((round(x$carapace.width[ix] / bin) * bin)[is.na(y)])
            r <- rbind(r, data.frame(x = as.numeric(names(tmp)), y = NA, n = as.numeric(tmp)))

            # Unconditional maturity identification:
            r$z <- rep(NA, nrow(r))
            r$z[which(r$x < 35)] <- 0 # Crab smaller than 30mm CW are considered immature.

            # Fit morphometric model:
            theta <- fit.morphometry.scsbio(r$x, r$y, r$z, n = r$n, sex = i, discrete = years[j] < 1998, trace = 0)

            # Extract maturity proportions:
            r$p <- morphometry.scsbio(r$x, r$y, r$z, theta = theta, fit = FALSE, discrete = years[j] < 1998)$p_mature_posterior

            # Map probabilities to original data:
            v[ix] <- r$p[gulf.utils::match(data.frame(x = round(x$carapace.width[ix] / bin) * bin, y = round(y / bin) * bin), r)]
         }
      }
   }

   # Convert to binary maturity:
   if (!probability) v <- (v >= 0.5)

   return(v)
}
