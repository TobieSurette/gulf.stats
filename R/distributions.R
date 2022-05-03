#' @title Probability Distributions
#'
#' @description Set of probability distributions.
#'
#' @param x Domain values.
#' @param n Number of random variates to generate.
#' @param eta Location parameter.
#' @param lambda Scale parameter.
#' @param gamma Skewness parameter parameter.
#' @param delta Shape parameter.
#' @param log Logical value specifying whether to return the logarithmic values of the function.
#'
#' @describeIn distributions Probability density function for the Johnson's SU distribution.
#' @export
djohnson <- function(x, eta = 0, lambda = 1, delta = 1, gamma = 0, log = FALSE){
   u <- (x - eta) / lambda
   v <- log(delta) - log(lambda) - 0.5 * log(2*pi) - 0.5 * log(1 + u*u) - 0.5 * (gamma + delta * asinh(u))^2

   if (!log) v <- exp(v)

   return(v)
}

#' @describeIn distributions Cumulative density function for the Johnson's SU distribution.
#' @export
pjohnson <- function(x, eta = 0, lambda = 1, delta = 1, gamma = 0, log = FALSE){
   u <- (x - eta) / lambda
   return(stats::pnorm(gamma + delta * asinh(u), log = log))
}

#' @describeIn distributions Cumulative density function for the Johnson's SU distribution.
#' @export
rjohnson <- function(n, eta = 0, lambda = 1, delta = 1, gamma = 0){
   return(lambda * sinh((qnorm(runif(n)) - gamma) / delta) + eta)
}
