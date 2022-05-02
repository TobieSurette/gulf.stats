#' @title Initial Touchdown and Lift-Off Times
#'
#' @description Functions to find initial parameter values for Minilog data.
#'
#' @param x Minilog data object.
#' @param truncate Touch
#'
#' @examples
#' x <- read.minilog(year = 2014, tow.id = tow.id, location = "tilt")

#' @rawNamespace S3method(fit,minilog)
fit.minilog <- function(x, truncate = FALSE){
   # Define log-likelihood function to find initial touchdown and liftoff points:
   loglike.mix <- function(theta, x, y, precision){
      # Regression model:
      mu    <- polynomial(theta[grep("beta", names(theta))])(x) # Regression mean.
      sigma <- exp(theta[["log.sigma"]])                        # Regression error for bottom phase.

      # Mixture density:
      p0 <- 1/ (1 + exp(theta[["logit.p"]]))
      d0 <- 1 / diff(range(y))                                                            # Uniform density.
      d1 <- pnorm(y + precision / 2, mu, sigma) - pnorm(y - precision / 2, mu, sigma) # Discrete Gaussian.
      d <-  p0 * d0 + (1-p0) * d1

      # Log-likelihood:
      v <- -sum(log(d))

      return(v)
   }

   # Define piecewise-linear likelihood:
   loglike.plm <- function(theta, x, y, fixed, precision){
      if (!missing(fixed)) theta <- c(theta, fixed)

      # Parse parameters:
      xp    <- theta[grep("^xp", names(theta))]   # Transition points.
      beta  <- theta[grep("beta", names(theta))]  # Bottom coeffcients.
      alpha_descent <- theta[grep("alpha_descent", names(theta))] # Descent phase basis function coefficients.
      alpha_ascent <- theta[grep("alpha_ascent", names(theta))]   # Ascent phase basis function coefficients.
      sigma <- exp(theta[grep("log.sigma", names(theta))])        # Regression error.

      # Define basis function for ascent and descent phases:
      basis <- function(x, beta){
         v <- polynomial(beta)(x)
         v[x < 0] <- 0
         return(v)
      }

      # Regression mean:
      mu <- basis(xp[1]-x, alpha_descent) + polynomial(beta)(x) + basis(x-xp[2], alpha_ascent)

      # Define outlier proportions:
      p <- as.numeric(1 / (1+exp(-theta["logit.p"])))
      p <- rep(p, length(x))
      p[(x >= xp[1]) & (x <= xp[2])] <- 1

      # Discrete likelihood:
      d <- p * (pnorm(y + precision / 2, mu, sigma[1]) - pnorm(y - precision / 2, mu, sigma[1])) +
           (1-p) * (pnorm(y + precision / 2, mu, sigma[1] + sigma[2]) - pnorm(y - precision / 2, mu, sigma[1] + sigma[2]))

      return(-sum(log(d)))
   }

   # Remove data tails:
   if (truncate){
      ix <- which(x$depth > 25)
      ix <- ix[round(length(ix) / 25):round(24*length(ix) / 25)]
      x <- x[ix, ]
   }
   if (nrow(x) == 0) return(as.POSIXct(c(NA, NA)))

   # Convert to seconds:
   x$t <- time2sec(time(x))
   x <- x[!is.na(x$depth) & !is.na(x$t), ]
   if (nrow(x) == 0) return(as.POSIXct(c(NA, NA)))

   # Initialize model parameters:
   tab <- table(round(x$depth,1))
   theta <- c(logit.p = 0, beta = as.numeric(names(tab)[which.max(tab)]), log.sigma = log(1))

   # Remove odd peaks:
   model <- mgcv::gam(x$depth ~ s(x$t))
   ix <- which(abs(residuals(model)) <= 4)
   #if (length(ix) > 0) x <- x[-ix, ]

   # Define measurement precision:
   precision <- table(abs(round(diff(x$depth),1)))
   precision <- precision[names(precision) != 0]
   precision <- as.numeric(names(precision)[which.max(precision)])

   # Fit constant model:
   for (i in 1:15) theta <- optim(theta, loglike.mix, x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(trace = 0))$par

   # Fit linear model:
   theta <- c(theta[-grep("beta", names(theta))], beta = as.numeric(c(theta[grep("beta", names(theta))], 0)))
   for (i in 1:15) theta <- optim(theta, loglike.mix,  x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(trace = 0))$par

   # Fit quadratic model
   #theta <- c(theta[-grep("beta", names(theta))], beta = as.numeric(c(theta[grep("beta", names(theta))], 0)))
   #for (i in 1:15) theta <- optim(theta, loglike.mix, x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(trace = 0))$par

   # Extract model parameters:
   xp <- theta[grep("^xp", names(theta))]
   beta <- theta[grep("beta", names(theta))]

   # Initial estimates for touchdown and liftoff:
   ix <- ix[which(x$depth[ix] >= (gulf.stats::polynomial(beta)(x$t[ix]) - 1))]

   # Initialize PLM parameters:
   theta <- c(xp = x$t[range(ix)],
              beta = as.numeric(coef(lm(x$depth[ix] ~ x$t[ix]))),
              alpha = c(-0.005, -0.005),
              log.sigma = c(3, 3),
              logit.p = 0)


   # Redefine data range:
   ix <- min(which((x$t < theta["xp1"]) & (x$depth > (polynomial(beta)(theta["xp1"]) - 40)))):max(which((x$t > theta["xp2"]) & x$depth > (polynomial(beta)(theta["xp2"]) - 40)))
   ix <- intersect(ix, which(x$depth > 25))

   fixed <- theta[c(grep("^xp", names(theta)), grep("beta", names(theta)))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   loglike.plm(theta, x = x$t, y = x$depth, precision = precision)
   loglike.plm(theta, x = x$t, y = x$depth, precision = precision, fixed = fixed)

   for (i in 1:5) theta <- optim(theta, loglike.plm, x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(trace = 3), fixed = fixed)$par

   theta <- c(fixed, theta)
   fixed <- theta[grep("^xp", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   for (i in 1:5) theta <- optim(theta, loglike.plm, x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(trace = 3), fixed = fixed)$par

   theta <- c(fixed, theta)
   for (i in 1:15) theta <- optim(theta, loglike.plm, x = x$t[ix], y = x$depth[ix], precision = precision, control =  list(maxit = 5000, trace = 3))$par


     # Parse parameters:
      xp    <- theta[grep("^xp", names(theta))]   # Transition points.
      beta  <- theta[grep("beta", names(theta))]  # Bottom coeffcients.
      alpha <- theta[grep("alpha", names(theta))] # Descent and ascent coefficients.

      # Define basis function for ascent and descent phases:
      basis <- function(x, degree = 2){
         v <- x^degree
         v[x < 0] <- 0
         return(v)
      }

      # Regression mean:
      mu <- alpha[1] * basis(xp[1]-x$t) + polynomial(beta)(x$t) + alpha[2] * basis(x$t-xp[2])

    #  mu <- -0.005 * basis(xp[1]-x$t) + polynomial(beta)(x$t) + -0.005 * basis(x$t-xp[2])

      plot(x$t, x$depth)
      vline(xp, col = "red", lwd = 2)
      lines(x$t, mu, lwd = 2, col = "blue")

   xp <- time(x[c(which.min(abs(x$t - xp[1])), which.min(abs(x$t - xp[2]))), ])

   return(xp)
}
