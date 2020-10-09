variogram.default <- function(x, z, lag, max.distance, mean.lag = TRUE, distance.units, response.units, trim, fit = TRUE, ...){
   # VARIOGRAM.DEFAULT - Fit a variogram object to data.
   #   'x' can be coordinates or a distance matrix.
   #   'z' : Vector of observations.
     
   # Check input arguments:
   if (missing(x)) stop ("Coordinate vector or matrix must be specified.")
   if (missing(z)) stop ("Observation vector must be specified.")
   if (!is.numeric(x)) stop ("'x' coordinate vector or matrix must be numeric.")
   if (!is.numeric(z)) stop ("'z' observation vector must be numeric.")
   if ((length(x) %% length(z)) != 0) stop("'x' and 'z' have incompatible lengths.")
   if (any(is.na(x) | is.na(z))) stop("'x' or 'z' contain NA values.")

   # Trim data:
   sigma <- sd(z)
   if (!missing(trim)){
      if (!is.null(trim) & !any(is.na(trim))){
         if (length(trim) == 1) if (trim > 0) trim <- -order(z)[(length(z)-trim+1):length(z)]
         if (length(trim) > 1) if (any(trim > 0)) stop("'trim' is not propoerly specified.")
         sigma <- sd(z)
         if (nrow(x) == ncol(x)) x <- x[trim, trim] else x <- x[trim, ]
         z <- z[trim]
      }
   }
   
   # Calculate distance matrix:
   if (is.vector(x)) x <- as.matrix(x)
   if ((nrow(x) == ncol(x)) & (all(x >= 0))) d <- as.matrix(x) else d <- as.matrix(dist(x))
   
   # Define maximum distance:
   if (missing(max.distance)) max.distance <- max(d) / 2
   if ((length(max.distance) != 1) | !is.numeric(max.distance)) stop("'max.distance' must be a numeric scalar.")
   if (max.distance > max(d)) stop(paste0("'max.distance' exceeds maximum observed distance between samples : ", max(d)))
   
   # Calculate distance matrix:
   index <- lower.tri(d) & d <= max.distance
   dd <- d[index]
   zz <- 0.5 * as.matrix(dist(z)) ^ 2 # Semi-variance.
   zz <- zz[index]

   # Calculate empirical variogram:
   if (!missing(lag)){
      # Calculate lag boundaries:
      lags <- seq(0, max.distance, by = lag)
      
      # Indices for lag intervals:
      ii <- floor(dd / lag) + 1 
      
      # Calculate summary results (empirical variogram):
      res  <- aggregate(list(lag = ii), list(lag = ii), unique)
      res <- data.frame(start.distance = lags[res$lag], end.distance = lags[res$lag+1] )
      if (!mean.lag){
         res$h <- (res$start.distance + res$end.distance) / 2
      }else{
         res$h <- aggregate(list(h = dd), list(lag = ii), mean)$h  # Calculate mean lag distance.
      }
      res$semi.variance <- aggregate(list(semi.variance = zz), list(lag = ii), mean)$semi.variance
      res$n <- aggregate(list(n = dd), list(lag = ii), length)$n # Number of pairs.
   }else{
      res <- data.frame(h = dd, semi.variance = zz, n = 1)
   }
   
   # Define 'lag':
   if (missing(lag)) lag <- NULL
   
   # Define object properties:
   v <- list(empirical = res,
             sd = sigma, 
             var = sigma^2,
             lag = lag,
             max.distance = max.distance)        
   
   if (!missing(distance.units)) v$distance.units <- distance.units
                
   # Define object class:
   class(v) <- c("variogram", "list")
   
   if (fit) v <- fit(v, ...)
   
   return(v)
}

