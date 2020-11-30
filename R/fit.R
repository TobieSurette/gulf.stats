#' Fit Statistical Models
#'
#' @description Functions to fit statistical models to data.
#'
#' @param x Data object or numeric vector representing a morphometric predictor variable..
#' @param y Numeric vector representing a morphometric response variable.
#' @param z Binary classification vector classifying morphometric data into discrete maturity groups.
#' @param n Number of observations for each (x,y) pair (optional).
#' @param sex Biological sex which specifies the type of analysis or initial values to be applied.
#' @param theta Parameter vector.
#' @param discrete Logical value specifying whether to treat observations as discrete rounded values.
#' @param print Logical value specifying whether to print information about progress during model fitting.
#' @param trace Optimization control output (see \code{\link[stats]{optim}}).
#' @param model Character string specifying the model type.
#' @param nugget Logical value specifying whether a variogram model contains a nugget semi-variance parameter.
#' @param distance.exponent Numeric value specifying the exponent to be applied in the distance metric.
#'

#' export
fit <- function(x, ...) UseMethod("fit")

# @describeIn fit Generic morphometric model fit method.
#' @export fit.morphometry
fit.morphometry <- function(x, ...) UseMethod("fit.morphometry")

#' @describeIn fit Fit a model to empirical variogram data.
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

#' @describeIn fit Fit a morphometric model to snow crab morphometric data.
#' @export fit.morphometry.scsbio
#' @rawNamespace S3method(fit.morphometry,scsbio)
fit.morphometry.scsbio <- function(x, y, z, n = 1, sex, theta, discrete = FALSE, print = FALSE, trace = 0){
   if (!missing(sex) & missing(theta)){
      if (sex == 1){
         # Male morphometry initial values:
         theta <- c(beta_immature = c(-2.03, 1.116, -0.06026, 0.0114), # Log-scale immature morphometric coefficients.
                    beta_mature = c(-2.858, 1.312),  # Log-scale mature morphometric coefficients.
                    log_sigma = -3.3,                # Log-scale standard error.
                    log_sigma_kurtotic = 0,          # Log-scale extra standard error for kurtotic observations.
                    logit_p_kurtotic = -2,           # Logit-scale proportion of kurtotic observations.
                    p_alpha = -11,                   # Logit-scale splm intercept parameter for mature proportions.
                    p_beta = c(0.25, 0.015, 0.25),   # Logit-scale splm slope parameters for mature proportions.
                    p_transition = c(45, 95),        # Logit-scale transition parameters for mature proportions.
                    p_window = 2.0)                  # Logit-scale splm window width parameter(s) for mature proportions.
      }
      if (sex == 2){
         # Female morphometry initial values:
         theta <- c(beta_immature = c(-2.72, 1.228), # Log-scale immature morphometric coefficients.
                    beta_mature = c(-2.80, 1.30),    # Log-scale mature morphometric coefficients.
                    log_sigma = -3,                  # Log-scale standard error.
                    log_sigma_kurtotic = 2,          # Log-scale extra standard error for kurtotic observations.
                    logit_p_kurtotic = -5.7,         # Logit-scale proportion of kurtotic observations.
                    p_alpha = -10.4,                 # Logit-scale splm intercept parameter for mature proportions.
                    p_beta = c(0.16, 0.015, 0.29),   # Logit-scale splm slope parameters for mature proportions.
                    p_transition = c(58.8, 101.1),   # Logit-scale transition parameters for mature proportions.
                    p_window = 1.45)                 # Logit-scale splm window width parameter(s) for mature proportions.
      }
   }

   # Negative log-likelihood function:
   loglike <- function(theta, x, y, z, n = 1, fixed, discrete = FALSE){
      if (!missing(fixed)) theta <- c(theta, fixed)
      v <- -sum(n * morphometry.scsbio(x, y, z, theta = theta, discrete = discrete)$loglike)
      return(v)
   }

   # Define optimization controls:
   control <- list(trace = trace, maxit = 1000)

   # Fit proportions
   if (print) cat("Fitting mature proportion parameters.\n")
   fixed <- theta[-grep("^p_", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, n = n, fixed = fixed, discrete = discrete, control = control)$par
   theta <- c(theta, fixed)

   # Fit kurtotic parameters:
   if (print) cat("Fitting kurtosis parameters.\n")
   fixed <- theta[-grep("kurtotic", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, n = n, fixed = fixed, discrete = discrete, control = control)$par
   theta <- c(theta, fixed)

   # Fit immature regression:
   if (print) cat("Fitting immature regression coefficients.\n")
   fixed <- theta[-grep("immature", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, n = n, fixed = fixed, discrete = discrete, control = control)$par
   theta <- c(theta, fixed)

   # Fit non-regression coefficients:
   if (print) cat("Fitting non-regression coefficients.\n")
   fixed <- theta[grep("mature", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, n = n, fixed = fixed, discrete = discrete, control = control)$par
   theta <- c(theta, fixed)

   # Fit immature regression:
   if (print) cat("Fitting complete model.\n")
   theta <- optim(theta, loglike, x = x, y = y, z = z, n = n, discrete = discrete, control = control)$par

   return(theta)
}
