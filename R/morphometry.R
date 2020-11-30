#' Analyze Morphometric Data
#'
#' @description Functions to analyze morphometric data.
#'
#' @param x Predictor variable.
#' @param y Response variable.
#' @param z Binary classification variable (optional).
#' @param theta Named model parameter vector.
#' @param sex Biological sex.
#' @param fit Logical value specifying whether to fit the model to the data.
#' @param discrete Logical value specifying whether the morphometric data are rounded.
#'
#' @examples
#' x <- read.scsbio(2010)
#' v <- morphometry(x)

#' @export
morphometry <- function(x, ...) UseMethod("morphometry")

#' @describeIn morphometry Perform morphometric regression analysis.
#' @export
morphometry.default <- function(x, y, species, kurtotic = FALSE, ...){
   # General case:
   m <- lm(log(y) ~ log(x))
   theta <- coef(m)
   names(theta) <- c("alpha", "beta")
   theta["alpha"] <- exp(theta["alpha"])

   return(theta)
}

#' @describeIn morphometry Perform morphometric regression analysis for snow crab biological data.
#' @export
#' @export morphometry.scsbio
morphometry.scsbio <- function(x, y, z, theta, sex, fit = TRUE, discrete = FALSE){
   if ("scsbio" %in% names(x)){
      x <- x$carapace.width
      y <- y$chela.height
      y[which(x$sex == 2)] <- y$abdomen.width
      sex <- x$sex
   }

   # Horner's method for polynomial evaluation::
   polyval <- function(x, p){
      p <- rev(p) # Intercept is last parameter:
      y <- p[[1]] * rep(1, length(x))
      if (length(p) > 1) for (i in 2:length(p)) y <- y * x + p[i]
      return(y)
   }

   # Calculate morphometric means:
   beta_immature <- theta[sort(names(theta)[grep("beta_immature", names(theta))])] # Immature coefficients.
   beta_mature   <- theta[sort(names(theta)[grep("beta_mature", names(theta))])]   # Mature coefficients.
   mu_immature   <- polyval(log(x), beta_immature)                                 # Immature mean.
   mu_mature     <- polyval(log(x), beta_mature)                                   # Mature mean.

   # Calculate unconditional maturity proportions:
   eta <- theta[grep("^p_", names(theta))] # Maturity proportion parameters.
   logit_p_mature <- splm(x, eta)
   p_mature <- 1 / (1 + exp(-logit_p_mature))

   # Error parameters:
   sigma <- exp(theta[["log_sigma"]])
   sigma_kurtotic <- exp(theta[["log_sigma_kurtotic"]])

   # Proportion of kurtotic distributional component:
   p_kurtotic <- 1 / (1 + exp(-theta[["logit_p_kurtotic"]]))

   # Compile results:
   v <- data.frame(mu_immature    = mu_immature,
                   mu_mature      = mu_mature,
                   p_mature       = p_mature,
                   sigma          = sigma,
                   sigma_kurtotic = sigma_kurtotic,
                   p_kurtotic     = p_kurtotic)

   # Calculate log-likelihood and posterior classification probabilities:
   if (!missing(y)){
       # Mixture component densities:
       if (discrete){
          dimm <- (1-p_kurtotic) * dnorm(log(y), mu_immature, sigma) + p_kurtotic * dnorm(log(y), mu_immature, sigma + sigma_kurtotic)
          dmat <- (1-p_kurtotic) * dnorm(log(y), mu_mature,   sigma) + p_kurtotic * dnorm(log(y), mu_mature,   sigma + sigma_kurtotic)
       }else{
          dimm <- ((1-p_kurtotic) * pnorm(log(y+0.5), mu_immature, sigma) + p_kurtotic * pnorm(log(y+0.5), mu_immature, sigma + sigma_kurtotic)) -
                  ((1-p_kurtotic) * pnorm(log(y-0.5), mu_immature, sigma) + p_kurtotic * pnorm(log(y-0.5), mu_immature, sigma + sigma_kurtotic))
          dmat <- ((1-p_kurtotic) * pnorm(log(y+0.5), mu_mature,   sigma) + p_kurtotic * pnorm(log(y+0.5), mu_mature,   sigma + sigma_kurtotic)) -
                  ((1-p_kurtotic) * pnorm(log(y-0.5), mu_mature,   sigma) + p_kurtotic * pnorm(log(y-0.5), mu_mature,   sigma + sigma_kurtotic))
       }

       # Initialize log-likelihood:
       v$loglike <- rep(0, length(x))

       # Likelihood for known classification:
       if (!missing(z)){
          ix <-!is.na(z)
          v$loglike[ix] <- (1-z[ix]) * log(1-p_mature[ix]) + z[ix] * log(p_mature[ix]) # Bernoulli likelihood.
          iy <- which(ix & (z == 0) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dimm[iy]) # Gaussian immature.
          iy <- which(ix & (z == 1) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dmat[iy]) # Gaussian mature.
       }
       ix <- is.na(x) | is.na(y)      # Missing morphometry.
       v$loglike[!ix]  <- v$loglike[!ix] + log((1-p_mature[!ix]) * dimm[!ix] + p_mature[!ix] * dmat[!ix])  # Add mixture contribution.

       # Posterior maturity probabilities:
       v$p_mature_posterior <- p_mature * dmat / ((1-p_mature) * dimm + p_mature * dmat)
       v$p_mature_posterior[is.na(v$p_mature_posterior)] <- p_mature[is.na(v$p_mature_posterior)]
   }

   return(v)
}
