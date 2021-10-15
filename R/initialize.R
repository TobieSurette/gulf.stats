#' @title Initialize Model Parameters
#'
#' @description Methods to provide initial parameter estimates for a given statistical model.
#'
#' @param x,y Data vectors.
#'
#' Initialize PLM parameters

#' @export
initialize <- function(x, ...) UseMethod("initialize")

#' @export
init <- initialize

#' @export
initialize.plm <- function(x, y, ...){
   # Mixture density function:
   dmix <- function(x, theta){
      mu <- theta[["mu"]]
      sigma <- exp(theta[["log_sigma"]])
      p <- 1 / (1 + exp(theta[["logit_p"]]))
      return(cbind((1-p) / (max(x) - min(x)), p * dnorm(x, mu, sigma)))
   }

   # Mixture density function:
   loglike.mix <- function(theta, x){
      v <- dmix(x, theta)
      v <- log(apply(v, 1, sum))
      return(-sum(v))
   }

   # Fit mixture model:
   theta <- c(mu = max(y), log_sigma = log(sd(y)/10), logit_p = 0)
   theta <- optim(theta, loglike.mix, x = y,  control = list(maxit = 10000, trace = 3))$par

   # Identify bottom points:
   v <- dmix(y, theta)
   ix <- (v[,2] / apply(v, 1, sum)) > 0.5

   # Derive initial values for bottom points:
   model <- lm(y[ix] ~ x[ix])
   theta <- c(slope2 = coef(model)[[2]],
              knots1 = x[ix][1], knots2 = x[ix][sum(ix)],
              log_sigma = log(sd(y[ix])))

   # Calculate touchdown linear parameters:
   ix <- x < theta[["knots1"]]
   model <- lm(y[ix] ~ x[ix])
   theta <- c(theta,
              intercept = coef(model)[[1]],
              slope1 = coef(model)[[2]])

   # Calculate liftoff linear parameters:
   ix <- x > theta[["knots2"]]
   model <- lm(y[ix] ~ x[ix])
   theta <- c(theta,
              slope3 = coef(model)[[2]])

   # Set default scale parameters:
   theta <- c(theta,
              log_scale1 = log(50),
              log_scale2 = log(50))

   # Order parameters:
   theta <- theta[sort(names(theta))]

   # Intercept correction:
   theta["intercept"] <- theta["intercept"] + mean(y - plm(theta = theta)(x))

   # Calculate PLM outer bounds:
   bounds <- c(theta[["knots1"]] - exp(theta[["log_scale1"]])/2,
               theta[["knots2"]] + exp(theta[["log_scale2"]])/2)

   theta["log_sigma_descent"] <- 0
   theta["log_sigma_ascent"]  <- 0
   theta["log_scale_descent"] <- log(40)
   theta["log_scale_ascent"]  <- log(40)

   # Define descent and ascent error parameters:
   ix <- x <= bounds[1]
   theta["log_sigma_descent"] <- log(sd(y[ix] - plm(theta = theta)(x[ix])))
   ix <- x >= bounds[2]
   theta["log_sigma_ascent"]  <- log(sd(y[ix] - plm(theta = theta)(x[ix])))

   return(theta)
}

