#' Kriging Spatial Modelling
#'
#' @description Spatial interpolation and integration funcions using kriging.
#'
#' @param x  Data coordinates or a data frame containing coordinates.
#' @param x0 Prediction coordinates.
#' @param y  Data response values.
#' @param z  External drift values at data locations.
#' @param z0 External drift values at prediction locations.
#' @param nugget Numerical value specifying the nugget value for the variogram model.
#' @param sill Numerical value specifying the sill or asymptotic value for the variogram model.
#' @param range Numerical value specifying the range value for the variogram model.
#' @param positive Logical value specifying whether to constrain kriged interpolation to be positive values.
#' @param variables Character string specifying the varigram in \code{x} to be analyzed.
#' @param variogram \code{variogram} object.
#' @param year Snow crab survey year.
#' @param category Snow crab survey year.
#' @param weight = FALSE,
#' @param hard.shelled = FALSE,
#' @param lag Lag distance, in kilometers, used for calculating the empirical variogram (i.e. bin size).
#' @param max.distance Maximum distance used for calculating the empirical variogram.
#' @param variogram.average Number of years to use in the variogram averaging.
#' @param polygon Coordinates of
#' @param grid Two-element vector specifying the dimensions of the interpolation prediction grid.
#' @param long.ref,lat.ref Lower-left coordinates used, along with \code{grid}, to define the interpolation grid over which kriging will be performed.
#' @param n Number of nearest neighbours used when kriging
#' @param nugget,sill,range Variogram model parameter values.
#' @param xlim,ylim Horizontal and vertical limits in kilometers with respect to the reference coordinates (\code{long.ref} and \code{lat.ref}).
#' @param polygon Polygon coordinates for within which kriging is to be performed.
#' @param bug Logical value specifying whether to include a bug in the calculations which was present in past analyses.

#' @export
ked <- function(x) UseMethod("ked")

#' @describeIn ked Perform kriging with external drift.
#' @export
ked.default <- function(x, x0, y, z, z0, nugget = 0, sill, range, positive = TRUE){
   # KED - Perform point kriging with external drift.
   #   'x'  : Data coordinates.
   #   'x0' : Prediction coordinates.
   #   'y'  : Data response values.
   #   'z'  : External drift values at data locations.
   #   'z0' : External drift values at prediction locations.

   # Number of data points:
   n <- nrow(x)

   # Basic functions:
   step <- function(x) return((x >= 0) + 1 - 1)
   spherical <- function(h) return((1-step(h-1)) * (1 - (1.5*h - .5*(h^3))))

   # Nugget value contribution to covariance matrix:
   k <- cbind(diag(rep(nugget, n)), 0)

   # Sill contribution to covariance matrix:
   t <- as.matrix(rbind(x, x0) / range)
   t <- t %*% t(t)
   h = sqrt(-2*t+diag(t) %*% matrix(1,ncol = n+1) + matrix(1, nrow = n+1) %*% t(diag(t)))
   h = h[1:n,]
   k <- k + sill * spherical(h)

   # Separate covariance matrix components:
   k0 = k[,n+1] # Prediction versus data covariance values.
   k = k[,1:n]  # Data covariance values.

   # ExtetrAdd external drift variable:
   k <- rbind(cbind(k, 1), 1)
   k[n+1, n+1] <- 0
   k <- rbind(cbind(k, c(z, 0)), c(z, 0, 0))
   k0 <- matrix(c(k0, 1), ncol = 1)
   k0 <- rbind(k0, z0)

   # Solve cokriging system by Gaussian elimination:
   const <- max(diag(k)) / max(max(abs(k[,(n+1):ncol(k)])))
   k[1:n,1:n] <- k[1:n,1:n] / const
   k0[1:n,] <- k0[1:n,] / const

   # Get kriging weights:
   w <- solve(k) %*% k0

   # calculation of cokriging estimates
   x0s = sum(w[1:n,] * y)

   # calculation of cokriging variances
   s = (nugget + sill) / const
   t = t(w) %*% k0
   s = const * (s-t);

   # Set negative values to zero:
   if (positive) x0s[x0s < 0] <- 0

   return(list(mean = x0s, var = s, weight = w, K = k, K0 = k0))
}

#' @describeIn ked Perform kriging with external drift for snow crab data.
#' @export
ked.scsset <- function(x, y, variables, variogram, year, category, weight = FALSE, hard.shelled = FALSE,
                       lag = 3, max.distance = 75, variogram.average = 1, xlim = c(0, 440), ylim = c(0, 400),
                       grid = c(100, 100), long.ref = -66, lat.ref = 45.5, n = 8, nugget, sill, range,
                       polygon = kriging.polygons()$gulf, bug = FALSE, ...){

   # Library for interpolating depth data.
   library(akima)

   # Append kriging response variables:
   if (!missing(y)){
      y <- as.data.frame(y)
      variables <- names(y)
      if (missing(x)) stop("'x' must be specified.")
      x <- cbind(x, y)
   }

   # Extract variogram parameters:
   parameters <- NULL
   if (!missing(sill) & !missing(range)){
      if (length(sill) != length(range)) stop("Variogram parameters must have the same number of elements.")
      if (missing(nugget)) nugget <- rep(0, length(sill))
      if (length(sill) != length(nugget)) stop("Variogram parameters must have the same number of elements.")
      parameters <- data.frame(nugget = nugget, sill = sill, range = range)
   }

   # Load appropriate data and variable and calculate variogram(s):
   if (!missing(year)){
      cat(paste0("Loading data for ", year, ".\n"))

      # Define variables
      if (missing(category) & !missing(variables)) category <- variables
      if (!is.null(parameters) | !missing(variogram)) variogram.average <- 1

      # Define number of years of data to be loaded.
      years <- (year - variogram.average + 1):year

      # Read survey data:
      x <- read.scsset(year = years, valid = 1)
      x <- x[(x$month >= 7), ]
      x <- x[substr(x$tow.id,2,2) != "C", ]
      if (missing(category)) stop("'category' must be specified.")

      # Calculate frequency observations:
      x <- summary(x, category = category, weight = weight, hard.shelled = hard.shelled, ...)

      # Scale to square kilometers:
      x[category] <- 1000000 * x[category] / repvec(x$swept.area, ncol = length(category))
      variables <- category
   }
   if (missing(variables) & !missing(category)) variables <- category

   years <- sort(unique(x$year))

   # Estimate variograms:
   if (is.null(parameters) & missing(variogram)){
      cat(paste0("Fitting variograms. \n"))
      if (variogram.average > 1){
         if (!('year' %in% names(x))) stop("'year' must be contained in the data to perform variogram averaging.")

         # Define matrix of fitted variograms:
         v <- matrix(list(), nrow = length(years), ncol = length(variables))

         for (j in 1:length(variables)){
            for (i in 1:length(years)){
               v[[i,j]] <- gulf::variogram(x[x$year == years[i], ], variable = variables[j], lag = 3, max.distance = max.distance, fit = TRUE, inits = list(range = 20))
            }
         }

         # Perform variogram averaging:
         w <- list()
         for (j in 1:length(variables)){
            res <- data.frame(start.distance = sort(unique(as.numeric(unlist(lapply(v[,j], function(x) x$empirical$start.distance))))))
            res$end.distance <- res$start.distance + lag
            for (k in 1:nrow(v)){
               tmp <- data.frame(name = rep(NA, nrow(res)))
               names(tmp) <- years[k]
               res <- cbind(res, tmp)
               index <- match(v[[k,j]]$empirical$start.distance, res$start.distance)
               res[index, ncol(res)] <- v[[k,j]]$empirical$semi.variance / v[[k,j]]$var
            }
            res$h <- NA
            res$n <- NA
            index <- match(v[[nrow(v),j]]$empirical$start.distance, res$start.distance)
            res[index, "h"] <- v[[nrow(v),j]]$empirical$h
            res[index, "n"] <- v[[nrow(v),j]]$empirical$n
            res$semi.variance <- v[[nrow(v),j]]$var * apply(res[, as.character(years)], 1, mean, na.rm = TRUE)
            if (length(which(is.na(res$h) | is.na(res$n))) > 0) res <- res[-which(is.na(res$h) | is.na(res$n)), ]

            w[[j]] <- v[[nrow(v),j]]
            w[[j]]$empirical$semi.variance <- res$semi.variance
            w[[j]] <- fit(w[[j]], distance.exponent = 1)
         }
         variogram <- w

         # Keep only recent data, not the data used for the averaging:
         x <- x[x$year >= min(x$year) + variogram.average - 1, ]
      }else{
         variogram <- list()
         print(variables)
         for (j in 1:length(variables)) variogram[[j]] <- variogram.scset(x, variable = variables[j], lag = 3, max.distance = max.distance, fit = TRUE, inits = list(range = 20))
         names(variogram) <- years
      }
   }

   # Extract variogram parameters:
   if (!missing(variogram)){
      if ("variogram" %in% class(variogram)) stop("'variogram' must be a 'variogram' object.")
      parameters <- as.data.frame(matrix(NA, nrow = length(variogram), ncol = 3))
      for (i in 1:length(variogram)){
         parameters[i,1] <- variogram[[i]]$nugget
         parameters[i,2] <- variogram[[i]]$sill
         parameters[i,3] <- variogram[[i]]$range
      }
      names(parameters) <- c("nugget", "sill", "range")
      rownames(parameters) <- variables
   }

   # Convert coordinates to kilometers:
   tmp <- deg2km(longitude(x), latitude(x), long.ref = lonref, lat.ref = latref, method = "ellipsoid")
   names(tmp) <- c("xkm", "ykm")
   x <- cbind(x, tmp)

   # Load kriging depth file:
   data(kriging.depth)

   # Interpolate depth grid at data locations:
   x$depth <- NA
   for (i in 1:nrow(x)){
      dx <- x$xkm[i] - depth$xkm
      dy <- x$ykm[i] - depth$ykm
      index <- order(dx^2 + dy ^2)[1:4]
      xx <- depth$xkm[index]
      yy <- depth$ykm[index]
      zz <- depth$depth[index]
      x$depth[i] <- interp(x = xx, y = yy, z = zz, x$xkm[i], x$ykm[i], linear = TRUE, extrap = TRUE, duplicate = "mean")$z
   }

   # Define kriging interpolation grid:
   if (all(grid == c(100, 100)) & (lonref == -66) & (latref == 45.5) & all(xlim == c(0, 440)) & all(ylim == c(0, 400))){
      data("kriging.grid100x100")
      x0 <- kriging.grid100x100
   }else{
      xx <- seq(xlim[1], xlim[2], by = diff(xlim)/(grid[1]-1))
      yy <- seq(ylim[1], ylim[2], by = diff(ylim)/(grid[2]-1))
      x0 <- expand.grid(xx, yy)
      names(x0) <- c("xkm", "ykm")
      x0 <- as.data.frame(x0)

      # Interpolate depth grid at kriging grid locations:
      x0$depth <- NA
      cat("Interpolating depth at kriging interpolation locations:\n")
      for (i in 1:nrow(x0)){
         if ((i %% 500) == 0) cat(paste0("   Interpolating ", i, " of ", nrow(x0), " points.\n"))
         dx <- x0$xkm[i] - depth$xkm
         dy <- x0$ykm[i] - depth$ykm
         index <- head(order(dx^2 + dy^2))[1:4]
         xx <- depth$xkm[index]
         yy <- depth$ykm[index]
         zz <- depth$depth[index]
         flag <- all(x0$xkm[i] > xx) | all(x0$xkm[i] < xx) | all(x0$ykm[i] > yy) | all(x0$ykm[i] < yy)
         if (!flag){
            if (all(zz == 0)) x0$depth[i] <- 0 else x0$depth[i] <- interp(x = xx, y = yy, z = zz, x0$xkm[i], x0$ykm[i], linear = TRUE, extrap = TRUE, duplicate = "mean")$z[1,1]
         }
      }

      # Append longitude and latitude coordinates:
      x0 <- cbind(x0, as.data.frame(km2deg(x0$xkm, x0$ykm, long.ref = -66, lat.ref = 45.5, method = "ellipsoid")))
   }

   near <- function(x, x0, n){
      # NEAR - Get indices of nearest data n points in each of the four cartesian quadrants.

      d = (x[,1]-x0[1])^2 + (x[,2]-x0[2])^2;

      ii <- list()
      ii[[1]] = which(x[,1] >= x0[1] & x[,2] >= x0[2])
      ii[[2]] = which(x[,1] < x0[1] & x[,2] >= x0[2])
      ii[[3]] = which(x[,1] < x0[1] & x[,2] < x0[2])
      ii[[4]] = which(x[,1] >= x0[1] & x[,2] < x0[2])

      for (i in 1:4){
         di = (x[ii[[i]],1]-x0[1])^2 + (x[ii[[i]],2]-x0[2])^2;
         ii[[i]] = ii[[i]][order(di)[1:min(n, length(di))]]
      }

      ii <- unlist(ii)
      ii <- ii[!is.na(ii)]

      return(ii)
   }

   # Recreate Matlab bug:
   if (bug){
      tmp <- parameters
      for (i in 1:nrow(parameters)){
         tmp$nugget[i] <- as.matrix(parameters)[1]
         tmp$sill[i] <- as.matrix(parameters)[2] + as.matrix(parameters)[1]
      }
      parameters <- tmp
   }

   # Perform kriging:
   cat("Performing kriging:\n")
   mu <- array(NA, dim = c(grid, length(variables)))
   var <- array(NA, dim = c(grid, length(variables)))
   for (i in 1:length(variables)){
      mu.tmp <- matrix(NA, nrow = nrow(x0))
      var.tmp <- mu.tmp
      cat(paste0("   Kriging '", variables[i], "' variable:\n"))
      for (j in 1:nrow(x0)){
         if ((j %% 500) == 0) cat(paste0("      Kriging ", j, " of ", nrow(x0), " points.\n"))

         id <- near(x[, c("xkm", "ykm")], unlist(x0[j, c("xkm", "ykm")]), n)

         res <- ked.default(x = x[id, c("xkm", "ykm")],
                            x0 = x0[j, c("xkm", "ykm")],
                            y = x[id, variables[i]],
                            z = x[id, "depth"],
                            z0 = x0[j, "depth"],
                            nugget = parameters$nugget[i],
                            sill   = parameters$sill[i] - parameters$nugget[i],
                            range  = parameters$range[i])

         mu.tmp[j] <- res$mean
         var.tmp[j] <- res$var
      }
      cat("\n")

      dim(mu.tmp) <- grid
      dim(var.tmp) <- grid
      mu[,,i] <- mu.tmp
      var[,,i] <- var.tmp
   }

   # Perform cross-validation:
   cv <- matrix(NA, nrow = nrow(x), ncol = length(variables))
   colnames(cv) <- variables
   for (i in 1:length(variables)){
      for (j in 1:nrow(x)){
         # Remove observation:
         xx <- x[-j, ]

         id <- near(xx[, c("xkm", "ykm")], unlist(x[j, c("xkm", "ykm")]), n)

         res <- ked.default(x = xx[id, c("xkm", "ykm")],
                            x0 = x[j, c("xkm", "ykm")],
                            y = xx[id, variables[i]],
                            z = xx[id, "depth"],
                            z0 = x[j, "depth"],
                            nugget = parameters$nugget[i],
                            sill   = parameters$sill[i] - parameters$nugget[i],
                            range  = parameters$range[i])

         cv[j,i] <- res$mean
      }
      cat("\n")
   }

   # Get kriging grid coordinates:
   mlon <- x0$longitude
   mlat <- x0$latitude
   dim(mlon) <- grid
   dim(mlat) <- grid

   # Compile results:
   res <- list(data = x,
               cross.validation = cv,
               variables = variables,
               variogram = variogram,
               map = mu,
               map.sigma = sqrt(var),
               map.longitude = mlon,
               map.latitude = mlat,
               dim = grid,
               reference = c(lonref, latref),
               xlim = xlim,
               ylim = ylim)

   class(res) <- c("ked", "list")

   return(res)
}
