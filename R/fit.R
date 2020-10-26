#' Fit Model to Data
#'
#' @description Functions to fit statistical models to data.

#' export
fit <- function(x, ...) UseMethod("fit")

#' @describeIn variogram Fit a model to empirical variogram data.
#' @export fit.variogram
#' @export
fit.variogram <- function(x, model = "spherical", nugget = TRUE, distance.exponent = 0, inits, ...){
   # Parse input arguments:
   model <- match.arg(tolower(model), c("exponential", "spherical", "gaussian"))

   # Define various variogram models:
   if (model == "exponential"){
      vfun <- function(h, nugget = 0, range = 1, sill = 1){
         v <- rep(sill, length(h))
         dim(v) <- dim(h)
         v <- (sill - nugget) * (1 - exp(-(3*h)/range)) + nugget
         return(v)
      }
   }
   if (model == "gaussian"){
      vfun <- function(h, nugget = 0, range = 1, sill = 1){
         v <- rep(sill, length(h))
         dim(v) <- dim(h)
         v <- (sill - nugget) * (1 - exp(-(3*(h^2))/(range^2))) + nugget
         return(v)
      }
   }
   if (model == "spherical"){
      vfun <- function(h, nugget = 0, range = 1, sill = 1){
         v <- rep(sill, length(h))
         dim(v) <- dim(h)
         index <- h < range
         v[index] <- (sill - nugget) * (((3 * h[index]) / (2* range)) - (h[index] ^ 3) / (2 * (range ^ 3))) + nugget
         return(v)
      }
   }

   # Extract empirical variogram values:
   res <- x$empirical

   # Find initial values for variogram parameters:
   if (missing(inits)) inits <- NULL
   if ("range" %in% names(inits)) range <- inits$range else range <- x$max.distance / 2
   names(range) <- "range"
   sill <- mean(res$semi.variance[round(nrow(res)/2):nrow(res)])
   names(sill) <- "sill"
   if (nugget){
      if (!is.null(x$lag)){
         index <- ((res$start.distance + res$end.distance) / 2) < (0.75 * range)
         tmp <- ((res$start.distance[index] + res$end.distance[index]) / 2)
         nugget <- coef(lm(res$semi.variance[index] ~ tmp, weights = res$n[index]))[1]
      }else{
         index <- res$h < (0.75 * range)
         nugget <- coef(lm(res$semi.variance[index] ~ res$h[index]))[1]
      }
      nugget <- max(c(0.0000001, nugget))
      names(nugget) <- "nugget"
   }else{
      nugget <- NULL
   }

   # Switch nugget and sill parameters if they are inverted:
   if (!is.null(nugget)){
      if (sill < nugget){
         tmp <- sill
         sill <- nugget
         nugget <- tmp
      }
      sill <- abs(sill - nugget)
      nugget <- abs(nugget)
   }else{
      sill <- abs(sill)
   }
   range <- abs(range)

   # Catenate parameter vector:
   theta <- c(nugget, range, sill)

   # Define objective function:
   ss <- function(theta, lag, semi.variance, weight, distance.exponent = 0){
      if ("nugget" %in% names(theta)) nugget <- abs(theta["nugget"]) else nugget <- 0
      sill <- nugget + abs(theta["sill"])
      range <- abs(theta["range"])
      mu <- vfun(lag, nugget = nugget, range = range, sill = sill)
      if (missing(weight)) weight <- 1
      if (distance.exponent != 0) weight <- weight / (lag ^ distance.exponent) # Add lag distance-weighted component.
      v <- sum(weight * (semi.variance - mu)^2)
      return(v)
   }

   # Estimate parameters:
   p <- optim(theta, ss, lag = res$h, semi.variance = res$semi.variance, weight = res$n , distance.exponent = distance.exponent, control = list(maxit = 4000, trace = 0))
   for (i in 1:5) p <- optim(p$par, ss, lag = res$h, semi.variance = res$semi.variance, weight = res$n , distance.exponent = distance.exponent, control = list(maxit = 4000, trace = 0))
   theta <- p$par

   # Parse parameter vector:
   if ("nugget" %in% names(theta)) nugget <- abs(theta["nugget"]) else nugget <- 0
   sill <- nugget + abs(theta["sill"])
   range <- abs(theta["range"])

   # Add parameters and model to object:
   x$nugget <- nugget
   x$sill   <- sill
   x$range  <- range
   x$model  <- model
   x$vfun   <- vfun

   return(x)
}
