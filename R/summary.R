#' Summary for Statistical Models
#'
#' @name summary
#'
#' @description Generate a summary for a statistical model.

#' @export
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
         mu <- x$map[,,j]
         ii <- (res$variable == x$variables[j]) & (res$polygon == names(polygon[i]))

         # Arithmetic statistics:
         id <- in.polygon(as.polygon(p$longitude, p$latitude), lon(scsset(x$data)), lat(scsset(x$data)))
         res$mean.sample[ii] <- mean(x$data[id,x$variables[j]]) * res$area[ii]
         res$sd.sample[ii] <- sd(x$data[id,x$variables[j]]) * res$area[ii]
         res$n.sample[ii] <- length(x$data[id,x$variables[j]])

         # Cross-validation stats:
         mad <- function(x) mean(abs(x - mean(x)))
         res$mean.CV[ii] <- mean(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])
         res$mad.CV[ii] <- mad(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])
         res$sd.CV[ii] <- sd(x$data[id,x$variables[j]] - x$cross.validation[id,x$variables[j]])

         # Global mean:
         res$mean[ii] <- mean(mu[index], na.rm = TRUE) * res$area[ii]

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

         id <- in.polygon(as.polygon(p$x, p$y), xx0[, 1], xx0[, 2])
         xx0 <- xx0[id,]

         # Data points which are close to the polygon boundary:
         tmp <- prcomp(cbind(p$x, p$y))
         u <- -tmp$rot
         l <- diag(tmp$sdev^2)
         co <- -tmp$x
         xc <- c(mean(p$x), mean(p$y))
         mi <- apply(co, 2, min)
         ma <- apply(co, 2, max)
         dd <- 0.3*sqrt((ma[1]-mi[1])*(ma[2]-mi[2]))
         cot <- rbind(mi-dd, c(mi[1]-dd, ma[2]+dd), ma+dd, c(ma[1]+dd, mi[2]-dd))
         poly2 <- cot %*% u + cbind(rep(xc[1], nrow(cot)), rep(xc[2], nrow(cot)))
         id <- in.polygon(as.polygon(poly2[, 1],poly2[, 2]), x$data$xkm, x$data$ykm)
         if (sum(id) == 0) id <- in.polygon(as.polygon(p$x, p$y), x$data$xkm, x$data$ykm)
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
            slog = sqrt(log((sigma^2)/(res$mean[ii]^2)+1));
            mlog = log(res$mean[ii])-(slog^2)/2;
            cilog = exp(c(mlog - 1.959964 * slog, mlog + 1.959964 *slog));
            res$lowerCI[ii] <- cilog[1]
            res$upperCI[ii] <- cilog[2]
         }
      }
   }

   return(res)
}
