#' Smoothed Piecewise Linear Models
#'
#' @description Functions for defining and evaluating smoothed piecewise-linear model.
#'
#' @param x Numeric vector at which the function is to be evaluated.
#' @param theta Named parameter vector with an intercept parameter \code{alpha}, slope parameters \code{beta0, ..., betak},
#'              where \code{k} is the number of transition (i.e. break) points, \code{transition1, ..., transitionk} and
#'              \code{window}, the transition window scale parameter (log-scale). \code{window} may also vary by transition
#'              point.
#'
#' @return If \code{x} is specified, a vector of function evaluations is returned, other wise the function itself is returned.
#'
#' @examples
#' # Single transition:
#' theta <- c(alpha = 0, slope = c(0, 1), transition = 0, window = -3)
#' x <- seq(-1, 1, length = 1000)
#' plot(x, splm(theta = theta)(x), type = "l", lwd = 2)
#' plot(x, splm(theta = theta, weights = TRUE)(x), type = "l", lwd = 2)
#'
#' # Multiple transitions:
#' theta <- c(alpha = 0, slope = c(0, 1, -1, 0), transition = c(0, 1, 2), window = -3)
#' x <- seq(-1, 3, length = 1000)
#' plot(x, splm(theta = theta)(x), type = "l", lwd = 2)

#' @export
splm <- function(x, ...) UseMethod("splm")

#' @describeIn splm Define or evaluate a smoothed piecewise-linear model.
#' @export
splm.default <- function(x, theta, weights = FALSE, ...){
   # Define function components:

   names(theta) <- tolower(names(theta))

   # Alias substitutions:
   names(theta) <- gsub("intercept", "alpha", names(theta))
   names(theta) <- gsub("slopes", "slope", names(theta))
   names(theta) <- gsub("slope", "beta", names(theta))
   names(theta) <- gsub("xp", "transition", names(theta))

   # Parse 'theta':
   alpha <- theta[sort(names(theta)[grep("alpha", names(theta))])]
   beta <- theta[sort(names(theta)[grep("beta", names(theta))])]
   transition <- theta[sort(names(theta)[grep("transition", names(theta))])]
   window <- exp(theta[sort(names(theta)[grep("window", names(theta))])])
   k <- length(transition) # Model order.
   if (length(window) == 1) window <- rep(window, k)

   # Define corresponding logistic weights for each component:
   if (weights){
      weight <- function(x){
         if (k == 0) return(rep(1,length(x)))
         eta <- matrix(NA, nrow = length(x), ncol = k)
         for (i in 1:k) eta[,i] <- (x - transition[[i]]) / window[[i]]
         if (k == 1) return(1 / (1 + exp(-eta)))
         return(NULL)
      }
   }

   # Define 'splm' function:
   mu <- function(x){
      y <- alpha + beta[1] * x
      for (i in 1:k) y <- y + window[i] * (beta[i+1] - beta[i]) * log(1 + exp((x - transition[i]) / window[i]))
      names(y) <- names(x)
      return(y)
   }

   # Return results:
   if (missing(x)){
      # Return functions:
      if (!weights) return(mu) else return(weight)
   }else{
      # Evaluate functions:
      if (!weights) return(mu(x)) else return(weight(x))
   }
}
