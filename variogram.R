#' Empirical Variograms
#'
#' @description Functions to
#'
#' @param x Data observations, coordinates or a distance matrix.
#' @param h Lag distances.
#' @param z Vector of data observations.
#' @param model Character string specifying a variogram function type.
#' @param nugget Logical value specifying whether the variogram model includes a nugget (i.e. zero-lag) variance component.
#' @param distance.exponent Positive numeric value specifying the exponent of the distance metric.
#' @param inits Initial parameter values.
#'

#' @describeIn variogram Generic \code{variogram} method.
#' @export
variogram <- function(x, ...) UseMethod("variogram")

#' @describeIn variogram Calculate empirical variogram.
#' @export
variogram.default <- function(x, z, lag, max.distance, mean.lag = TRUE, distance.units, response.units, trim, fit = TRUE, ...){
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

#' @describeIn variogram Calculate empirical variogram for snow crab survey data.
#' @export
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

#' @describeIn variogram Evaluate variogram at specified distances.
#' @export
variogram.variogram <- function(x, h)  v <- x$vfun(h, nugget = x$nugget, sill = x$sill, range = x$range)

#' @describeIn variogram Fit a model to empirical variogram data.
#' @export fit.variogram
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

#' @describeIn variogram Plot a variogram object.
#' @export
plot.variogram <- function(x, scale, cex = 1.5, cex.axis = 1.25, xaxt = "s", yaxt = "s",
                           xlab, ylab, xlim, ylim, show.nugget = TRUE, show.sill = TRUE, show.range = TRUE, ...){
   # PLOT.VARIOGRAM - Graphically display a variogram.

   # Extract model parameters:
   res <- x$empirical

   # Determine y-scaling:
   v <- floor(mean(log10(res$semi.variance)))
   if (missing(scale)){
      if (v > 1) scale <- 10^v else scale <- 1
   }

   print(scale)
   # Define axis labels:
   if (missing(xlab)) xlab <- "Distance"
   if (!is.null(x$distance.units) & (xlab != "")) xlab <- paste0(xlab, " (", x$distance.units, ")")
   if (missing(ylab)) ylab <- "Semi-variance"
   if ((scale != 1) & (ylab != "")){
      a <- floor(mean(log10(scale)))
      ylab <- parse(text = paste("Semi-variance (10^{", a, "})"))
   }

   h <- seq(0, x$max.distance, len = 1000)

   if (missing(xlim)) xlim <- c(0, x$max.distance)
   if (missing(ylim)) ylim <- c(0, 1.1*max(res$semi.variance / scale))

   plot(res$h, res$semi.variance / scale, pch = 21, bg = "grey",
        xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxs = "i", yaxs = "i",
        cex = cex, cex.axis = cex.axis, type = "n", xaxt = "n", yaxt = "n", ...)
   mtext(xlab, 1, 2.5, cex = cex)
   mtext(ylab, 2, 2, cex = cex)
   grid()

   lines(par("usr")[1:2], c(x$var, x$var) / scale, col = "chartreuse3", lwd = 2, lty = "dashed")
   if (xaxt != "n") axis(1)
   if (yaxt != "n") axis(2)

   # Plot model and parameters:
   if (!is.null(x$sill)){
      nugget <- x$nugget
      sill <- x$sill
      range <- x$range
      lines(h, x$vfun(h, nugget = nugget, sill = sill, range = range) / scale, lwd = 2, col = "red")
      if (show.sill) lines(par("usr")[1:2], c(sill, sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (show.range) lines(c(range, range), c(0, sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (show.nugget) lines(par("usr")[1:2], c(nugget, nugget) / scale, lwd = 1, lty = "dashed", col = "red")
   }

   # Display empirical semi-variances:
   points(res$h, res$semi.variance / scale, pch = 21, bg = "grey", cex = 1.5)
   if (!is.null(x$lag)) text(res$h, res$semi.variance / scale, res$n, pos = 3, cex = 0.75)

   # Display parameter texts:
   #if (!is.null(x$sill)){
      #if (show.range){
      #   if (range > mean(par("usr")[1:2])) pos = 2 else pos = 4
      #   text(range, 0.5 * sill / scale, paste0("range = ", round(range,2), " ", x$distance.units), cex = 1.25, pos = pos)
      #}
      #if (show.sill) text(0, sill / scale + 0.025*diff(par("usr")[3:4]), paste0("sill = ", round((sill) / scale,3)), cex = 1.25, pos = 4)
      #if (show.nugget & (nugget / sill) > 0.01) text(0, nugget / scale - 0.025*diff(par("usr")[3:4]), paste0("nugget = ", round(nugget / scale,3)), cex = 1.25, pos = 4)
   #}
   #text(0, (x$var / scale) + 0.030*diff(par("usr")[3:4]), paste0("Var(sample) = ", round(x$var / scale,3)), cex = 1.25, pos = 4, col = "chartreuse3")
   #print(str(x))
   #v <- floor(log10(par("usr")[4]))
   str <- ""
   if (is.finite(range) & !is.na(range))
   str <- c(str, paste0("range = ", round(range,1), x$distance.units))
   str <- c(str, paste0("sill = ", round(sill/ 10^v ,3), " x 10^", v),
                 paste0("nugget = ", round(nugget / 10^v,3), " x 10^", v),
                 paste0("Var[z] = ", round(x$var / 10^v,3), " x 10^", v))

   #parse(text = )

   dx <- diff(par("usr")[1:2])
   dy <- diff(par("usr")[3:4])
   for (i in 1:length(str)){
      text(par("usr")[1] + 0.55 * dx,
           par("usr")[3] + 0.30 * dy - 0.06 * (i-1) * dy,
           str[i], pos = 4)
   }

   box()
}

