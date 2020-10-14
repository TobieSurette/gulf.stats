#' Kriging
#'
#' @description
#'
#'

ked.default <- function(x, x0, y, z, z0, nugget, sill, range, positive = TRUE){
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

ked.scset <- function(x, y, variables, variogram, year, category, weight = FALSE, hard.shelled = FALSE, lag = 3, max.distance = 75, variogram.average = 1,
                      xlim = c(0, 440), ylim = c(0, 400), grid = c(100, 100), lonref = -66, latref = 45.5, n = 8, nugget, sill, range,
                      polygon = kriging.polygons()$gulf, bug = FALSE, ...){
   # KED.SCSET - Perform kriging with external drift for snow crab data.

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
      x <- read.scset(year = years, valid = 1)
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

summary.ked <- function(x, polygon){
   # Generate a summary of a 'ked' object.

   # Parse polygon argument:
   if (is.numeric(polygon)){
      if (ncol(polygon) != 2) stop("'polygon' must have two columns.")
      if (is.null(colnames(polygon))) colnames(polygon) <- c("longitude", "latitude")
      polygon <- as.data.frame(polygon)
   }
   if (is.data.frame(polygon)) polygon <- list(polygon)
   for (i in 1:length(polygon)){
      polygon[[i]] <- as.data.frame(polygon[[i]])
      str <- tolower(colnames(polygon[[i]]))
      str[grep("lon", str)] <- "longitude"
      str[grep("lat", str)] <- "latitude"
      colnames(polygon[[i]]) <- str
   }

   # Convert coordinates to kilometers:
   areas <- NULL
   for (i in 1:length(polygon)){
      p <- polygon[[i]]
      p <- cbind(p[setdiff(names(p), c("x", "y"))], deg2km(p$longitude, p$latitude, long.ref = x$reference[1], lat.ref = x$reference[2], method = "ellipsoid"))
      polygon[[i]] <- p
      areas[i] <- area.polygon(as.polygon(p$x, p$y))
   }

   # Convert polygon coordinates to kilometers:
   res <- expand.grid(x$variables, 1:length(polygon))
   names(res) <- c("variable", "polygon")
   res <- as.data.frame(res)
   res[, 1] <- as.character(res[, 1])
   res$area <- areas[res$polygon]
   if (!is.null(names(polygon))) res$polygon <- names(polygon)[res$polygon]
   res$mean.sample <- NA
   res$sd.sample <- NA
   res$n.sample <- NA
   res$mean <- NA
   res$sd <- NA
   res$lowerCI <- NA
   res$upperCI <- NA
   res$mean.CV <- NA
   res$mad.CV <- NA
   res$sd.CV <- NA
   for (i in 1:length(polygon)){
      p <- polygon[[i]]
      for (j in 1:length(x$variables)){
         # Calculate global mean:
         index <- in.polygon(as.polygon(p$longitude, p$latitude), x$map.longitude, x$map.latitude)
         dim(index) <- x$dim
         mu <- x$map[,,j]
         ii <- (res$variable == x$variables[j]) & (res$polygon == names(polygon[i]))

         # Arithmetic statistics:
         id <- in.polygon(as.polygon(p$longitude, p$latitude), longitude.scset(x$data), latitude.scset(x$data))
         res$mean.sample[ii] <- mean(x$data[id,x$variables[j]]) *  res$area[ii]
         res$sd.sample[ii] <- sd(x$data[id,x$variables[j]]) *  res$area[ii]
         res$n.sample[ii] <- length(x$data[id,x$variables[j]])

         # Cross-validation stats:
         mad <- function(x) mean(abs(x - mean(x)))
         res$mean.CV[ii] <- mean(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])
         res$mad.CV[ii] <- mad(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])
         res$sd.CV[ii] <- sd(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])

         # Global mean:
         global <- mean(mu[index], na.rm = TRUE) * res$area[ii]
         res$mean[ii] <- global

         # Calculate global variance:
         area <- res$area[ii]
         ndisc <- 800
         lx <- diff(range(p$x))
         ly <- diff(range(p$y))
         nbdis <- floor(ndisc * lx * ly / area)+1;
         nx <- floor(sqrt(nbdis))+1;

         gx <- seq(min(p$x), max(p$x), by = lx / nx)
         gy <- seq(min(p$y), max(p$y), by = ly / nx)
         xx0 <- expand.grid(gx, gy)

         id <- in.polygon(as.polygon(p$x, p$y), xx0[, 1], xx0[, 2]); # inpolygon est une fonction de matlab 5.
         xx0 <- xx0[id,]

         # Data points which are close to the polygon boundary:
         tmp <- prcomp(cbind(p$x, p$y));
         u <- -tmp$rot
         l <- diag(tmp$sdev^2)
         co <- -tmp$x
         xc <- c(mean(p$x), mean(p$y));
         mi <- apply(co, 2, min);
         ma <- apply(co, 2, max);
         dd <- 0.3*sqrt((ma[1]-mi[1])*(ma[2]-mi[2]));
         cot <- rbind(mi-dd, c(mi[1]-dd, ma[2]+dd), ma+dd, c(ma[1]+dd, mi[2]-dd));
         poly2 <- cot %*% u + cbind(rep(xc[1], nrow(cot)), rep(xc[2], nrow(cot)))
         id <- in.polygon(as.polygon(poly2[, 1],poly2[, 2]), x$data$xkm, x$data$ykm);
         if (sum(id) == 0) id <- in.polygon(as.polygon(p$x, p$y), x$data$xkm, x$data$ykm);
         if (sum(id) == 0){
             xx <- repvec(p$x, ncol = length(x$data$xkm)) - repvec(x$data$xkm, nrow = length(p$x))
             yy <- repvec(p$y, ncol = length(x$data$ykm)) - repvec(x$data$ykm, nrow = length(p$y))
             dd <- sqrt(xx^2 + yy^2)
             id <- apply(dd, 1, function(x) order(x)[1:5])
             dim(id) <- NULL
             id <- sort(unique(id))
         }
         data <- x$data[id, c("xkm", "ykm", "depth")]

         # Define covariance model:
         step <- function(x) return((x >= 0) + 1 - 1)
         spherical <- function(h) return((1-step(h-1)) * (1 - (1.5*h - .5*(h^3))))

         # Distance matrices:
         hyy <- sqrt((repvec(xx0[,1], nrow = nrow(xx0)) - repvec(xx0[,1], ncol = nrow(xx0)))^2 +
                    (repvec(xx0[,2], nrow = nrow(xx0)) - repvec(xx0[,2], ncol = nrow(xx0)))^2)

         # Coviarance and cross-covariance matrices:
         Cyy <- (x$variogram[[j]]$sill - x$variogram[[j]]$nugget) * spherical(hyy / x$variogram[[j]]$range) # Prediction covariance.
         names(xx0) <- colnames(data[, 1:2])
         cx = rbind(data[,c("xkm", "ykm")], xx0)
         k <- matrix(0, nrow = nrow(data), ncol = nrow(cx))
         ccx <- cx / x$variogram[[j]]$range
         h <- sqrt((repvec(ccx$xkm, nrow = nrow(ccx)) - repvec(ccx$xkm, ncol = nrow(ccx)))^2 +
                   (repvec(ccx$ykm, nrow = nrow(ccx)) - repvec(ccx$ykm, ncol = nrow(ccx)))^2)
         h <- h[1:nrow(data), ]
         for (kk in 1:nrow(data)) k[kk,kk] <- x$variogram[[j]]$nugget
         k <- k + (x$variogram[[j]]$sill - x$variogram[[j]]$nugget) *  spherical(h)
         k0 = k[,(nrow(data)+1):ncol(k)]
         k = k[,1:nrow(data)]
         k <- cbind(rbind(k, 1), 1)
         k[nrow(k), nrow(k)] <- 0
         k0 <- rbind(k0, 1)
         k0 <- apply(k0, 1, mean)
         if (rcond(k) > 1e-15){
            l = solve(k) %*% k0;

            # Global standard error:
            sigma <- area * sqrt(mean(Cyy) - sum(l[, 1] * k0))
            res$sd[ii] <- sigma

            # Lognormal confidence intervals:
            slog = sqrt(log((sigma^2)/(global^2)+1));
            mlog = log(global)-(slog^2)/2;
            cilog = exp(c(mlog - 1.959964 * slog, mlog + 1.959964 *slog));
            res$lowerCI[ii] <- cilog[1]
            res$upperCI[ii] <- cilog[2]
         }
     }
   }

   return(res)
}
