#' @title Estimate Trawl Event Times
#'
#' @description Functions to find initial parameter values for Star Oddi data.
#'
#' @param x Star Oddi data object.
#' @param t Time vector.
#' @param theta Named free parameter vector.
#' @param fixed Named fixed parameter vector.
#' @param model Model to be used for modelling bottom trawling phase.
#' @param y Pressure observations.
#'s
#' @examples
#' x <- read.star.oddi(year = 2014, tow.id = tow.id, location = "tilt")

init.tilt <- function(x, init){
   s <- read.scsset(year = unique(year(x)), tow.id = tow.id(x))
   init <- c(time(s, "start"), time(s, "stop"))

   # Isolate data:
   ix <- which((time(x) >= (init[1] - 90)) & (time(x) <= (init[2] + 90)))
   x <- x[ix, ]

   # Estimate kurtotic normal mixture parameters:
   theta <- gulf.stats::fit.knorm(x$tilt.x)

   # Classification probabilities:
   u <- theta["p"] * dnorm(x$tilt.x, theta["mu"], theta["sigma_0"])
   v <- (1-theta["p"]) * dnorm(x$tilt.x, theta["mu"], theta["sigma_1"])
   p <- u / (u + v)

   # Extract times:
   xp <- time(x)[range(which(p > 0.5))]
   names(xp) <- c("touchdown", "liftoff")

   return(xp)
}

loglike.tilt <- function(theta, t, x, fixed){
   if (!missing(fixed)) theta <- c(theta, fixed)

   # Parse model parameters:
   xp_touchdown <- theta[["xp_touchdown"]]
   xp_liftoff <- theta[["xp_liftoff"]]
   window_touchdown <- exp(theta[["log_window_touchdown"]])
   window_liftoff <- exp(theta[["log_window_liftoff"]])
   mu <- theta[grep("^mu", names(theta))]
   sigma <- exp(theta[setdiff(grep("^log_sigma", names(theta)), grep("^log_sigma_outlier", names(theta)))])

   logit_p_outlier <- theta[grep("logit_p_outlier", names(theta))]

   # Name model parameters:
   names(mu) <- gsub("mu_", "", names(mu))
   names(sigma) <- gsub("log_sigma_", "", names(sigma))

   # Define transition functions:
   p_touchdown <- 1 / (1 + exp(-((t - xp_touchdown) / window_touchdown)))
   p_liftoff   <- 1 / (1 + exp(-((t - xp_liftoff) / window_liftoff)))

   # Initialize likelihood matrix:
   d <- matrix(NA, nrow = nrow(x), ncol = 3)

   # Define variables:
   vars <- c("x", "y", "z")
   phases <- c("descent", "bottom", "ascent")

   for (j in 1:length(vars)){
      # Define mixture density components:
      d_descent <- dnorm(x[,paste0("tilt.", vars[j])], mu[paste0("descent", "_", vars[j])], sigma[paste0("descent", "_", vars[j])])
      d_bottom  <- dnorm(x[,paste0("tilt.", vars[j])], mu[paste0("bottom", "_", vars[j])],  sigma[paste0("bottom", "_", vars[j])])
      d_ascent  <- dnorm(x[,paste0("tilt.", vars[j])], mu[paste0("ascent", "_", vars[j])],  sigma[paste0("ascent", "_", vars[j])])

      # Mixture density:
      d[,j] <- d_descent + p_touchdown * (d_bottom - d_descent) + p_liftoff * (d_ascent - d_bottom)
   }

   return(-sum(log(d)))
}

loglike.pressure <- function(theta, x, y){
   # Remove missing values:
   x <- x[!is.na(y)]
   y <- y[!is.na(y)]

   # Data range:
   ylim <- range(y)

   # Mixture parameters:
   p <- 1/ (1 + exp(theta[["logit.p"]]))      # Mixture proportions.

   # Regression model:
   beta  <- theta[grep("beta", names(theta))] # Bottom polynomial parameters.
   mu    <- polynomial(beta)(x)               # Regression mean.
   sigma <- exp(theta[["log.sigma"]])         # Regression error for bottom phase.

   # Mixture likelihood:
   d <- (1-p) / diff(ylim) + p * dnorm(y, mu, sigma)
   v <- -sum(log(d))

   return(v)
}

event.times.pressure <- function(x, model = "linear"){
   # Determine bottom model degree:
   degree <- match(tolower(model), c("constant", "linear", "quadratic", "cubic"))

   # Truncate data out:
   #x <- x[, ]
   x <- x[(x$pressure > 1) & (x$pressure > (0.7 * max(x$pressure))), ]

   # Initialize model parameters:
   tab <- table(round(x$pressure,1))
   theta <- c(logit.p = 0,
              beta = as.numeric(names(tab)[which.max(tab)]),
              log.sigma = log(1))

   # Convert to seconds:
   x$t <- time2sec(time(x))

   # Fit mixture model:
   for (j in 1:degree){
      for (i in 1:15) theta <- optim(theta, loglike.pressure, x = x$t, y = x$pressure, control =  list(trace = 0))$par
      theta <- c(theta[-grep("beta", names(theta))], beta = c(theta[grep("beta", names(theta))], 0))
   }

   # Extract model parameters:
   beta <- theta[grep("beta", names(theta))]
   sigma <- exp(theta[grep("log.sigma", names(theta))])

   # Initial estimates for touchdown and liftoff:
   ix <- which((x$pressure >= (gulf.stats::polynomial(beta)(x$t) - 1.96 * sigma) & (x$pressure <= (polynomial(beta)(x$t) + 1.96 * sigma))))
   xp <- time(x[range(ix), ])
   names(xp) <- c("touchdown", "liftoff")

   return(xp)
}

event.times.tilt <- function(x, plot = FALSE, add = FALSE){
   # Initialize change-point values:
   #xp <- init.tilt(x)
   xp <- event.times.pressure(x)

   # Truncate data:
   t <- time(x)
   x <- x[which(t >= (xp[1] - 180) & t <= (xp[2] + 180)), ]
   t <- time(x)

   # Initial parameter values:
   vars <- c("tilt.x", "tilt.y", "tilt.z")
   phases <- c("descent", "bottom", "ascent")
   bounds <- c(time(x[1,]), xp, time(x[nrow(x),]))
   mu <- matrix(NA, nrow = length(phases), ncol = length(vars))
   sigma <- mu
   for (i in 1:length(phases)){
      for (j in 1:length(vars)){
         ix <- which((t >= bounds[i]) & (t <= bounds[i+1]))

         mu[i,j] <- mean(x[ix, vars[j]])
         sigma[i,j] <- sd(x[ix, vars[j]])
      }
   }
   dimnames(mu) <- list(phase = phases, variable = vars)
   dimnames(sigma) <- list(phase = phases, variable = vars)

   # Define relative time:
   reftime <- time(x[1, ])
   xp <- time2sec(xp, reftime)
   t <- time2sec(time(x), reftime)

   # Build parameter vector:
   theta <- c(xp_touchdown = as.numeric(xp)[1],
              xp_liftoff = as.numeric(xp)[2],
              log_window_touchdown = 2,
              log_window_liftoff = 2,
              mu_descent_x = as.numeric(mu["descent", "tilt.x"]),
              mu_descent_y = as.numeric(mu["descent", "tilt.y"]),
              mu_descent_z = as.numeric(mu["descent", "tilt.z"]),
              mu_bottom_x = as.numeric(mu["bottom", "tilt.x"]),
              mu_bottom_y = as.numeric(mu["bottom", "tilt.y"]),
              mu_bottom_z = as.numeric(mu["bottom", "tilt.z"]),
              mu_ascent_x = as.numeric(mu["ascent", "tilt.x"]),
              mu_ascent_y = as.numeric(mu["ascent", "tilt.y"]),
              mu_ascent_z = as.numeric(mu["ascent", "tilt.z"]),
              log_sigma_descent_x = log(as.numeric(sigma["descent", "tilt.x"])),
              log_sigma_descent_y = log(as.numeric(sigma["descent", "tilt.y"])),
              log_sigma_descent_z = log(as.numeric(sigma["descent", "tilt.z"])),
              log_sigma_bottom_x = log(as.numeric(sigma["bottom", "tilt.x"])),
              log_sigma_bottom_y = log(as.numeric(sigma["bottom", "tilt.y"])),
              log_sigma_bottom_z = log(as.numeric(sigma["bottom", "tilt.z"])),
              log_sigma_ascent_x = log(as.numeric(sigma["ascent", "tilt.x"])),
              log_sigma_ascent_y = log(as.numeric(sigma["ascent", "tilt.y"])),
              log_sigma_ascent_z = log(as.numeric(sigma["ascent", "tilt.z"])),
              logit_p_outlier = 0)

   v.old <- loglike.tilt(theta, t = t, x = x)
   for (i in 1:30){
      cat(paste(i, round(v.old, 4), "\n"))
      v <- optim(theta, loglike.tilt, t = t, x = x, control = list(trace = 0, maxit = 1500))
      v.old <- v$value
      theta <- v$par
   }

   # Parse model parameters:
   xp_touchdown <- theta[["xp_touchdown"]]
   xp_liftoff <- theta[["xp_liftoff"]]
   window_touchdown <- exp(theta[["log_window_touchdown"]])
   window_liftoff <- exp(theta[["log_window_liftoff"]])
   mu <- theta[grep("^mu", names(theta))]
   sigma <- exp(theta[setdiff(grep("^log_sigma", names(theta)), grep("^log_sigma_outlier", names(theta)))])

   # Name model parameters:
   names(mu) <- gsub("mu_", "", names(mu))
   names(sigma) <- gsub("log_sigma_", "", names(sigma))

   # Define transition functions:
   logistic <- function(x){
      v <- 1 / (1 + exp(-x))
      v[x > 300] <- 1
      v[x < -300] <- 0
      return(v)
   }
   p_touchdown <- logistic((t - xp_touchdown) / window_touchdown)
   p_liftoff   <- logistic((t - xp_liftoff) / window_liftoff)

   # ======================================================================================
   if (plot){
      plot(t, x$tilt.x, pch = 21, bg = "grey", cex = 0.5, ylim = c(-15, 60), yaxs = "i", add = add)
      grid()
      vline(theta[grep("xp", names(theta))], lwd = 2, col = "red", lty = "dashed")

      m <- mu["descent_x"] +
         p_touchdown * (mu["bottom_x"] - mu["descent_x"]) +
         p_liftoff * (mu["ascent_x"] - mu["bottom_x"])

      s <- sigma["descent_x"] +
         p_touchdown * (sigma["bottom_x"] - sigma["descent_x"]) +
         p_liftoff * (sigma["ascent_x"] - sigma["bottom_x"])

      lines(t, m, col = "blue", lwd = 3)
      lines(t, m - 1.96 * s, col = "blue", lwd = 2, lty = "dashed")
      lines(t, m + 1.96 * s, col = "blue", lwd = 2, lty = "dashed")
      hline(0, col = "red", lwd = 3, lty = "solid")
   }

   xp <- time(x[c(which.min(abs(t - theta["xp_touchdown"])), which.min(abs(t - theta["xp_liftoff"]))), ])

   return(xp)
}




