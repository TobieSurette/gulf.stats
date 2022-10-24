#' @title  Calculate Biomass
#'
#' @description Functions to calculate snow crab biomass.
#'
#' @param year Survey year.
#' @param category Snow crab category.
#' @param weight Logical value specifying whether to convert catches to weights.
#'
#' @examples
#' r <- biomass.grid(2022)
#' r <- biomass.grid(2022, category = "MMSC12") # Fishery recruits.

#' @export biomass.grid
biomass.grid <- function(year, category = "COM"){
   survey <- gulf.spatial::read.gulf.spatial("scs bounds")
   survey <- gulf.graphics::as.polygon(survey[,1], survey[,2])

   # Read three years of data (for variogram averaging):
   s <- gulf.data::read.scsset(year = (year-2):year, survey = "regular", valid = 1) # Tow data.
   b <- gulf.data::read.scsbio(year = (year-2):year, survey = "regular")            # Biological data.
   b$tow.id <- gulf.data::tow.id(b)

   # Import catch data:
   gulf.utils::import(s, fill = 0) <- gulf.data::catch(b, category = category, weight = TRUE, as.hard.shelled = TRUE, units = "t") # Merge catches.
   s[category] <- 1000000 * s[category] / repvec(s$swept.area, ncol = length(category))   # Convert to tonnes per km2.

   # Perform kriging:
   m <- ked(s, variables = category, variogram.average = 3, lag = 3, max.distance = 75)

   # Interpolation grid:
   x <- seq(-66.5, -60, len = 1000)
   y <- seq(45, 49, len = 1000)

   # Interpolate:
   ix <- !is.na(m$map.longitude) & !is.na(m$map.latitude) & !is.na(m$map[,,1])
   zz <- akima::interp(x = m$map.longitude[ix],
                       y = m$map.latitude[ix],
                       z = m$map[,,1][ix],
                       xo = x, yo = y,
                       linear = TRUE, extrap = TRUE, duplicate = "mean")$z

   xx <- t(repvec(x, nrow = length(y)))
   yy <- t(repvec(y, ncol = length(x)))
   ix <- gulf.graphics::in.polygon(survey, xx, yy)
   dim(ix) <- dim(zz)
   zz[!ix] <- NA

  # clg()
#   image(x, y, zz, col = cols(length(breaks)-1), breaks = c(seq(0, 10, by = 1), 25))

#   # Find which 10x10 grids overlap the survey polygon:
#   lines(survey[[1]]$x, survey[[1]]$y, col = "red")

   x0 <- expand.grid(seq(-67, -59, by = 1/6), seq(44, 50, by = 1/6))
   names(x0) <- c("longitude", "latitude")
   x0 <- x0[gulf.graphics::in.polygon(survey, x0[,1], x0[,2]), ]
   r <- cbind(x0[,1], x0[,2], x0[,1]+1/6, x0[,2]+1/6)
   r <- rbind(r, cbind(x0[,1]-1/6, x0[,2], x0[,1], x0[,2]+1/6))
   r <- rbind(r, cbind(x0[,1]-1/6, x0[,2]-1/6, x0[,1], x0[,2]))
   r <- rbind(r, cbind(x0[,1], x0[,2]-1/6, x0[,1]+1/6, x0[,2]))
   r <- unique(r)

   labels <- NULL
   for (i in 1:nrow(r)){
      labels <- c(labels, deg2grid((r[i,1]+r[i,3])/2, (r[i,2] + r[i,4]) / 2))
   }
   r <- aggregate(r, by = list(grid = labels), mean)

   ix <- !is.na(zz)
   xxx = xx[ix]
   yyy = yy[ix]
   zzz = zz[ix]
   r$biomass <- NA
   for (j in 1:nrow(r)){
      gx <- c(r$V1[j], r$V1[j], r$V3[j], r$V3[j], r$V1[j])
      gy <- c(r$V2[j], r$V4[j], r$V4[j], r$V2[j], r$V2[j])
      p <- gulf.graphics::as.polygon(gx, gy)
      ix <- which(gulf.graphics::in.polygon(p, xxx, yyy))
      tmp <- gulf.spatial::deg2km(gx, gy)

      xp <- seq(r$V1[j], r$V3[j], len = 50)
      yp <- seq(r$V2[j], r$V4[j], len = 50)
      xp <- expand.grid(xp, yp)
      proportion <- sum(gulf.graphics::in.polygon(survey, xp[,1], xp[,2])) / nrow(xp)

      area <- proportion * gulf.graphics::area(gulf.graphics::as.polygon(tmp[,1], tmp[,2]))
      #print(area)
      r$biomass[j] <- mean(zzz[ix]) * area
   }

   r$category <- category
   r$year <- year

   return(r[c("year", "category", "grid", "biomass")])
}
