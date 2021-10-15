#' @title Model Log-likelihood functions.
#'
#' @description Methods for calculating log-likelihood values for various types of models.
#'
#' @param theta Parameter vector.
#' @param x,y Data vectors.
#' @param fixed Vector of fixed parameter values.
#' @param precision Precision of observed data values.

#' @export
loglike <- function(x, ...) UseMethod("loglike")

#' @export log.likelihood
log.likelihood <- loglike

#' @export
loglike.plm <- function(theta, x, y, fixed, precision){
   # Append fixed parameters:
   if (!missing(fixed)){
      theta <- theta[setdiff(names(theta), names(fixed))]
      theta <- c(theta, fixed)
   }

   # Predictor funtion:
   mu <- plm(theta = theta)

   # Calculate error scaling function:
   sigma <- sigma.plm(theta = theta)

   # Residuals:
   eps <- y - mu(x)

   # AR correlation parameter:
   rho <- 1 / (1 + exp(-theta[["logit_rho"]]))

   # Variance of the first observation:
   sigma_0 <- sqrt((sigma(x[1])^2  +  exp(theta[["log_sigma"]])^2) / (1 - rho^2))

   # Calculate likelihood:
   if (missing(precision)){
      v <- rep(NA, length(x))
      v[1] <- dnorm(eps[1], 0, sigma_0, log = TRUE)
      v[2:length(eps)] <- dnorm(eps[-1], rho * eps[-length(eps)], sigma(x[2:length(eps)]), log = TRUE)
   }else{
      v[1] <- pnorm(eps[1] + precision / 2, 0, sigma_0) - pnorm(eps[1] - precision / 2, 0, sigma_0)
      v[2:length(eps)] <- dnorm(eps[-1], rho * eps[-length(eps)], sigma(x[2]), log = TRUE)
     # v <- log(pnorm(delta + precision/2, 0, sigma(x), log = FALSE) - pnorm(delta - precision/2, 0, sigma(x), log = FALSE))
   }

   return(-sum(v))
}
