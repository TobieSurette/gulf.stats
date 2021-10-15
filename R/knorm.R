#' @title Kurtotic Normal Mixture
#'
#' @description Fits a kurtotic normal mixture model to data.
#'
#' @param x Data vector.
#'

#' @export knorm
knorm <- function(x, mu, sigma, p){
   v <- (1-p) * dnorm(x, mu, sigma[2]) + p * dnorm(x, mu, sigma[1])
   return(v)
}

#' @export fit.knorm
fit.knorm <- function(x){
   # Define mixture log-likehood function:
   loglike <- function(theta, x){
      p <- 1/ (1 + exp(-theta[["logit.p"]]))   # Mixture proportions.
      mu <- theta[["mu"]]
      sigma <- as.numeric(exp(theta[grep("log.sigma", names(theta))]))
      sigma[2] <- sigma[1] + sigma[2]
      v <- knorm(x, mu, sigma, p)
      return(-sum(log(v)))
   }

   # Initial parameter estimates:
   theta <- c(mu = median(x), log.sigma = c(log(sd(x) / 5), log(sd(x) / 2)), logit.p = 0)

   # Fit model:
   theta <- optim(theta, loglike, x = x, control = list(trace = 0))$par

   # Parse parameters:
   p <- 1/ (1 + exp(-theta[["logit.p"]]))   # Mixture proportions.
   mu <- theta[["mu"]]
   sigma <- as.numeric(exp(theta[grep("log.sigma", names(theta))]))
   sigma[2] <- sigma[1] + sigma[2]
   theta <- c(p = p, mu = mu, sigma_0 = sigma[1], sigma_1 = sigma[2])

   return(theta)
}
