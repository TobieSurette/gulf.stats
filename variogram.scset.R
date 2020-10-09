variogram.scset <- function(x, category, variable, weight = FALSE, lag = 3, max.distance = 75, distance.exponent = 1, ...){
   # VARIOGRAM.SCSET - Calculate an empirical variogram for an 'scset' object.
    
   if (!missing(category)){
      # Read biological data:
      b <- read.scbio(year = unique(x$year))
      vars <- c("year", "tow.id")
      res <- aggregate(is.category(b, category), by = x[c("year", "tow.id")], sum, na.rm = TRUE)
      index <- match(res[vars], x[vars])
      x$number.caught <- NA
      x$number.caught[index] <- res$x
      x$density <- x$number.caught / x$swept.area  
      z <- x$density    
   }
   
   if (!missing(variable)) z <- x[, variable]
                                                    
   # Calculate distance matrix:
   d <- as.matrix(dist(deg2km(abs(longitude(x)), latitude(x), long.ref = 66, lat.ref = 45.5, method = "ellipsoid")))
   
   # Fit variogram to data:
   v <- variogram(d, z, lag = lag, max.distance = max.distance, model = "spherical", distance.units = "km", distance.exponent = distance.exponent, ...)
      
   return(v)
}
