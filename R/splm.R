#' Smoothed Piecewise Linear Models
#'
#' @name splm
#'
#' @description Piecewise linear models with smoothed transitions.
#'
#' @param order Model order, i.e. number of piecewise linear segments minus one.
#' @param intercept Intercept parameter of the first linear segment.
#' @param slope Slope parameters of each linear segment.
#' @param transition Knot locations.
#' @param window Window width parameters at each knot location.
#' @param window.type Window width parameterization.
#' @param sigma Error parameters for each linear segment.
#' @param transition.basis Transition kernel basis.
#' @param type Type of piecewise segment.
#' @param sigma.model Error model type.
#' @param family Error distribution.

#' @rdname splm
#' @export splm
splm <- function(x, ...) UseMethod("splm")

#' @rdname splm
#' @export
splm.default <- function(x, y, order, intercept, slope, transition, window, window.type, sigma,
                 transition.basis = "logistic", type = "linear", error.model = "constant",
                 family = "Gaussian", ...){

   # SPLM - Smoothed piecewise linear model.

   # Parse input data:
   if (!missing(x)){
      if (!is.numeric(x) | !is.vector(x)) stop("'x' must be a numeric vector.")
   }
   if (!missing(y)){
      if (missing(x)) stop("'x' must be specified.")
      if (!is.numeric(y) | !is.vector(y)) stop("'y' must be a numeric vector.")
   }else{
      if (!missing(x)) stop("'y' must be specified.")
   }
   if (!missing(x)){
      if (length(x) != length(y)) stop("'x' and 'y' must be the same length.")
   }else{
      x <- NULL
      y <- NULL
   }

   # Parse input parameters:
   transition.basis <- match.arg(tolower(transition.basis), c("step", "linear", "logistic"))
   type <- match.arg(tolower(type), c("constant", "linear"))
   error.model <- match.arg(tolower(error.model), c("constant", "exponential", "expconst", "logisticconst", "splm"))
   if (!missing(window)){
      if (length(window) > 1) window.type <- "variable"
      if (length(window) == 1) window.type <- "constant"
   }
   if (!missing(window.type)){
      window.type <- match.arg(tolower(window.type), c("constant", "variable"))
   }else{
      window.type <- "constant"
   }

   # Parse 'order' argument:
   if (missing(order)){
      if (!missing(transition)) order = length(transition)
      if (missing(order) & (type == "constant")){
         if (!missing(order) & !missing(intercept)){
            if ((length(intercept)-1) != order)
               stop("Model order and length of 'intercept' parameters is inconsistent.")
         }
         if (missing(order) & !missing(slope)) order <- length(slope) - 1
      }
      if  (missing(order) & (type == "linear")){
         if (!missing(order) & !missing(slope)){
            if ((length(slope)-1) != order)
               stop("Model order and length of 'slope' parameters is inconsistent.")
         }
         if (missing(order) & !missing(slope)) order <- length(slope) - 1
      }
      if (!missing(order) & !missing(window)){
         if (!((length(window) == 1) | (length(window) == order)))
            stop("Model order and length of 'window' parameters is inconsistent.")
      }
      if (missing(order) & !missing(window)){
         if (length(window) > 1){
            order <- length(window)
         }
      }
   }
   if (missing(order)) order <- 1

   # Check that order is a non-negative integer:
   if ((length(order) != 1) | !is.numeric(order))
      stop("'order' must be a non-negative integer.")
   if (((order %% 1) != 0) | (order < 0)) stop("'order' must be a non-negative integer.")
   if (missing(order)) stop("'order' must be defined for properly initialize 'theta' parameter vector.")

   # Build 'splm' object:
   v <- list(x = x,
             y = y,
             theta = NULL,
             order = order,
             transition.basis = transition.basis,
             type = type,
             window.type = window.type,
             error.model = error.model)
   class(v) <- unique(c("splm", class(v)))

   # Build parameter vector:
   if (v$type == "constant") theta <- parameter(n = v$order+1, label = "intercept")
   if (v$type == "linear") theta <- c(parameter(n = 1, label = "intercept"), parameter(n = v$order+1, label = "slope"))
   theta <- c(theta, parameter(n = v$order, label = "transition", link = "ordered"))
   if (v$transition.basis != "step"){
      if (v$window.type == "constant")  theta <- c(theta, parameter(n = 1, label = "window", link = "log"))
      if (v$window.type == "variable") theta <- c(theta, parameter(n = v$order, label = "window", link = "log"))
   }
   if (v$error.model == "constant")      theta <- c(theta, parameter(n = 1, label = "sigma"))
   if (v$error.model == "exponential")   theta <- c(theta, parameter(n = 2, label = "sigma"))
   if (v$error.model == "expconst")      theta <- c(theta, parameter(n = 3, label = "sigma"))
   if (v$error.model == "logisticconst") theta <- c(theta, parameter(n = 4, label = "sigma"))
   if (v$error.model == "splm")          theta <- c(theta, parameter(n = x$order+1, label = "sigma"))

   # Consistency checks:
   if ((v$transition.basis == "step") & !missing(window))
      stop("'window' parameter may only be specified when transition basis is a step function.")
   if ((v$type == "constant") & !missing(slope))
      stop("'slope' parameters may not be specified when 'type' is 'constant'.")

   # Assign parameter components:
   if (!missing(intercept))  theta(theta, "intercept")                    <- intercept
   if (!missing(slope))      theta(theta, "slope")                        <- slope
   if (!missing(transition)) theta(theta, "transition", transform = TRUE) <- transition
   if (!missing(window))     theta(theta, "window", transform = TRUE)     <- window
   if (!missing(sigma))      theta(theta, "sigma")                        <- sigma

   v$theta <- theta

   return(v)
}

print.splm <- function(x){
   # PRINT.SPLM - Print 'splm' object to the R console.

   # Function to generate properly formatted strings:
   bracket <- function(v, quote = "'", sort = FALSE, digits = 4){
      if (sort & is.numeric(v)) if (!any(is.na(v))) v <- sort(v)
      if (is.numeric(v)) v <- signif(v, digits)
      if (length(v) > 10){
         n <- length(v)
         if (is.character(v)) v <- paste(quote, v[c(1, length(v))], quote, sep = "")
         v <- paste("[", v[1], ", ..., ", v[length(v)], "]", sep = "")
      }else{
         if (is.character(v)) v <- paste(quote, v, quote, sep = "")
         if (length(v) > 1) v <- paste("[", paste(v, collapse = ", "), "]", sep = "")
      }
      return(v)
   }

   # Display summary of 'alkey' object:
   cat(paste0("'splm' object:\n"))

   # Display parameter values:
   theta <- theta(x, full = TRUE)
   cat(paste0("                x : ", bracket(x$x), "\n"))
   cat(paste0("                y : ", bracket(x$y), "\n"))
   cat(paste0("            order : ", bracket(x$order), "\n"))
   cat(paste0(" transition basis : ", bracket(x$transition.basis), "\n"))
   cat(paste0("             type : ", bracket(x$type), "\n"))
   cat(paste0("      error model : ", bracket(x$error.model), "\n"))
   cat(paste0("        intercept : ", bracket(theta(x, "intercept")), "\n"))
   if (x$type != "constant") cat(paste0("            slope : ", bracket(theta(x, "slope")), "\n"))
   transition <- as.vector(theta(x, "transition", transform = TRUE) )
   cat(paste0("       transition : ", bracket(transition)), "\n")
   if (x$transition.basis != "step"){
      window <- theta(x, "window", transform = TRUE)
      cat(paste0("           window : ", bracket(window)), "\n")
   }
   cat(paste0("            sigma : ", bracket(theta(x, "sigma")), "\n"))
   if ("prior" %in% names(x)){
      cat(paste0("         prior(s) :"))
      for (i in 1:length(x$prior)){
         str <- paste0(rep(" ", 18-nchar(names(x$prior[[i]]))), collapse = "")
         str <- paste0(str, "'", names(x$prior[i]), "':")
         str <- paste0(str, x$prior[[i]]$str, "\n")
         cat(str)
      }
   }
   cat(paste0("            theta : ", bracket(theta(x)), "\n"))
}

summary.splm <- function(x){

   v <- list()
   v$order <- x$order
   v$intercept <- theta(x, "intercept")
   v$slope <- theta(x, "slope")
   v$transition <- theta(x, "transition", transform = TRUE)
   v$controls <- predict(x, v$transition, smooth = FALSE)

   return(v)
}

prior.splm <- function(x, index, full = FALSE, transform = FALSE){
   # PRIOR.SPLM - Extract priors from  'splm' object.

   # Define 'index' argument:
   if (missing(index) & full) index <- rep(TRUE, length(x$theta))
   if (missing(index) & !full){
      if (!is.null(x$fixed)) index <- !x$fixed else index <- rep(TRUE, length(x$theta))
   }

   # Index is a character vector of parameter names:
   if (is.character(index)){
      if (all(tolower(index) %in% names(x$theta))){
         index <- match(tolower(index), names(x$theta))
         if (length(index) > 1) stop("Invalid value for variable names.")
      }else{
         index <- grep(tolower(index), tolower(names(x$theta)))
      }
      if (length(index) == 0) stop("Unable to match parameter names.")
   }

   # Index is a numeric vector:
   if (is.numeric(index)){
      if (any(index %% 1 != 0) | any(index < 1) | any(index > length(x$theta)))
         stop("Invalid index vector for the parameter vector.")
      tmp <- rep(FALSE, length(x$theta))
      tmp[index] <- TRUE
      index <- tmp
   }

   # Parameter vector to be returned:
   v <- x$prior

   return(v[index])
}
"prior<-.splm" <- function(x, value, index, full = FALSE, transform = FALSE){
   # PRIOR<-.SPLM - Parameter assignement method for an 'splm' object.

   # Reset parameter vector:
   if (missing(index)){
      if (is.null(value)) x$prior <- list()
      if (length(value) == 1) if (is.na(value)) x$prior <- list()
   }

   # Index is a character vector of parameter names:
   if (is.character(index)){
      if (all(tolower(index) %in% names(x$theta))){
         index <- match(tolower(index), names(x$theta))
         if (length(index) != length(value))
            stop("index and length of parameter vector are inconsistent.")
         value <- value[order(index)]
         index <- index[order(index)]
      }else{
         if (length(index) > 1) stop("Invalid value for vairable names.")
         index <- grep(tolower(index), tolower(names(x$theta)))
      }
      if (length(value) == 0) stop("Unable to match parameter names.")
   }

   # Index is a numeric vector:
   if (is.numeric(index)){
      if (any(index %% 1 != 0) | any(index < 1) | any(index > length(x$theta)))
         stop("Invalid index vector for the parameter vector.")
      tmp <- rep(FALSE, length(x$theta))
      tmp[index] <- TRUE
      index <- tmp
   }

   # Check that value is logical:
   if (!is.logical(index)) stop("Invalid index parameter specification.")
   if (any(x$fixed & index)) stop("Some target assigned parameters are fixed parameters.")

   # Check that 'value' and 'index' are consistent:
   if (length(value) == 1) value <- rep(value, sum(index))
   if (!((length(value) == sum(index)) | (length(value) == length(index))))
      stop("index and length of parameter vector are inconsistent.")

   # Check theta value has the proper format:
   if (!is.character(value)) stop("'value' must be a character string.")
   #value <- "N(0,45)"
   str <- value
   str <- gsub("[(),]", " ", str)
   str <- gsub(" +", " ", str)
   str <- strsplit(str, " ")[[1]]
   if (!(str[1] %in% c("N", "LN"))) stop("Invalid prior distribution specification.")
   if (any(gsub("[0-9.]", "", str[2:length(str)]) != "")) stop("Prior parameters must be numeric.")

   # Define prior functions:
   tmp <- list(str = value)
   if (str[1] == "N")
      if (length(str) == 3)
         tmp$fun <- function(x) return(dnorm(x, as.numeric(str[2]), as.numeric(str[3])))
   if (str[1] == "LN")
      if (length(str) == 3)
         tmp$fun <- function(x) return(dlnorm(x, as.numeric(str[2]), as.numeric(str[3])))
   tmp <- list(tmp)
   names(tmp) <- names(x$theta)[index]

   print(tmp)

   # Initialize 'prior' field:
   if (is.null(x$prior)) x$prior <- list()

   # Assign value:
   if (!(names(tmp) %in% names(x$prior))) x$prior <- c(x$prior, tmp)

   return(x)
}

predict.splm <- function(x, x0, type = "mean", smooth = TRUE){
   # PREDICT.SPLM - Returns the predicted value associated with an 'splm'.

   # Parse 'x0' argument:
   if (missing(x0)) x0 <- x$x

   # Parse input arguments:
   type <- match.arg(tolower(type), c("mean", "mu", "sigma", "error"))
   if (type == "mu")    type <- "mean"
   if (type == "error") type <- "sigma"

   # Integrated step function:
   if (!smooth | (x$transition.basis == "step")){
      fun <- function(x){
         v <- rep(0, length(x))
         v[x > 0] <- x[x > 0]
         return(v)
      }
   }else{
      # Integrated linear (i.e. quadratic) transtion:
      if (x$transition.basis == "linear"){
         fun <- function(x){
            v <- 0.5 * (x * (x + 1) + 0.25)
            v[x < -0.5] <- 0
            v[x > 0.5] <- x[x > 0.5]
            return(v)
         }
      }

      # Integrated logistic transtion:
      if  (x$transition.basis == "logistic"){
         fun <- function(x){
            v <- 0.25 * log(1 + exp(4*x))
            v[x > 150] <- x[x > 150]
            return(v)
         }
      }
   }

   # Extract full parameter vector:
   theta <- theta(x, full = TRUE, transform = TRUE)
   intercept <- theta[grep("intercept", names(theta))]
   slope <- theta[grep("slope", names(theta))]
   transition <- theta[grep("transition", names(theta))]

   # Expand window parameter:
   if (x$transition.basis != "step"){
      window <- theta[grep("window", names(theta))]
      if (length(window) == 1) window <- rep(window, x$order)
   }

   # Calculate error prediction:
   if (type == "sigma"){
      sigma <- theta[grep("sigma", names(theta))]
      if (x$error.model == "constant")      v <- exp(sigma[1] * rep(1, length(x0)))
      if (x$error.model == "exponential")   v <- exp(sigma[1] +  sigma[2] * x0)
      if (x$error.model == "expconst")      v <- exp(sigma[1]) + exp(sigma[2] * (x0 - sigma[3]))
      if (x$error.model == "logisticconst") v <- exp(sigma[1]) + exp(sigma[4]) * (1 / (1 + exp(-sigma[2] -  sigma[3] * x0)))
      if (x$error.model == "splm"){
         v <- rep(sigma[1], length(x0))
         for (i in 1:x$order){
            if (x$transition.basis == "logistic"){
               v <- v + (sigma[i+1] - sigma[i]) * (1 / (1 + exp(-4*(x0 - transition[i]) / window[i])))
            }
         }
      }
      return(v)
   }

   # Piecewise constant model:
   if (x$type == "constant"){
      v <- rep(intercept[1], length(x0))
      for (i in 1:x$order){
         v <- v + window * (slope[i+1]-slope[i]) * fun((x0 - transition[i]) / window[i])
      }
   }

   # Define linear model for the first component:
   if (x$type == "constant") v <- rep(intercept[1], length(x0))
   if (x$type == "linear")   v <- intercept[1] + slope[1] * x0

   # Define other linear components:
   if (!smooth | x$transition.basis == "step"){
      # Piecewise constant model:
      for (i in 1:x$order){
         v <- v + (slope[i+1]-slope[i]) * fun((x0 - transition[i]))
      }
   }else{
      # Piecewise linear model:
      for (i in 1:x$order){
         v <- v + window[i] * (slope[i+1]-slope[i]) * fun((x0 - transition[i]) / window[i])
      }
   }

   # Strip names:
   v <- as.vector(v)

   return(v)
}

simulate.splm <- function(x, n = 1, x0, discrete = FALSE, drop = TRUE, ...){
   # SIMULATE.SPLM - Generate random samples from an 'splm' object.

   # Define predictor variable:
   if (missing(x0)) x0 <- x$x

   # Observation 'error':
   sigma <- predict(x, x0, sigma = TRUE, ...)

   # Regression mean:
   mu <- predict(x, x0, ...)

   # Generate random samples:
   v <- matrix(nrow = length(x0), ncol = n)
   for (i in 1:n)  v[,i] <- mu + sigma * rnorm(length(x0))

   # Round off simulated response data the same as the observed responses:
   if (discrete){
      p <- precision(s$y)
      v <- round(v / p) * p
   }

   # Convert result to vector form:
   if (drop & any(dim(v)==1)) v <- as.vector(v)

   return(v)
}

residuals.splm <- function(x){
   # RESIDUALS.SPLM - Calculate reaiduals for an 'splm' object.

   mu <- predict(x)
   sigma <- predict(x, sigma = TRUE)
   r <- (x$y - mu) / sigma

   return(r)
}

init.splm <- function(x){
   # INIT.SPLM - Initialize parameter vector of an 'splm' object.

   # Initialize transition parameters:
   xp <- theta(x, "transition", transform = TRUE)
   if (any(is.na(xp))){
      xp <- quantile(x$x, seq(0, 1, len = length(xp)+2))
      xp <- xp[2:(length(xp)-1)]
   }

   # Initialize linear coeffcients:
   b <- theta(x, "intercept")
   if (x$type == "constant"){
      if (any(is.na(b))){
         tmp <- c(min(x$x), xp, max(x$x))
         for (i in 1:(length(xp)+1)){
            index <- (x$x >= tmp[i]) & (x$x <= tmp[i+1])
            if (sum(index) > 1) b[i] <- mean(x$y[index])
         }
      }
   }else{
      m <- theta(x, "slope", transform = TRUE)
      if (any(is.na(c(b, m)))){
         tmp <- c(min(x$x), xp, max(x$x))
         for (i in 1:(length(xp)+1)){
            index <- (x$x >= tmp[i]) & (x$x <= tmp[i+1])
            if (sum(index) > 2){
               r <- lm(x$y[index] ~ x$x[index])
               if (i == 1) b <- coef(r)[1]
               m[i] <- coef(r)[2]
            }
         }
      }
   }

   # Initialize window parameter(s):
   if (x$transition.basis != "step"){
      w <- theta(x, "window", transform = TRUE)
      if (any(is.na(w))){
         w <- rep(0.1 * diff(range(x$x)) / (x$order+1), length(w))
      }
   }

   # Initialize 'sigma':
   sigma <- theta(x, "sigma")
   if (length(sigma) == 1){
      if (names(sigma) == "sigma") names(sigma) <- "sigma0"
   }
   if (is.na(sigma["sigma0"])){
      sigma["sigma0"] <- log(sd(x$y))
      sigma[setdiff(names(sigma), "sigma0")] <- 0
   }

   # Assign parameters:
   theta(x, "intercept") <- b
   if (x$type != "constant") theta(x, "slope") <- m
   theta(x, "transition", transform = TRUE) <- xp
   if (x$transition.basis != "step") theta(x, "window", transform = TRUE) <- w
   theta(x, "sigma") <- sigma

   return(x)
}

# THETA.SPLM - Extract 'splm' parameter vector.
theta.splm <- function(x, ...) return(theta(x$theta, ...))

loglike.splm <- function(theta, x, discrete = FALSE, ...){
   # LOG.LIKELIHOOD.SPLM - Log-likelihood function for an 'splm' object.

   # Parse input parameters:
   if ("splm" %in% class(theta)){
       x <- theta
   }else{
      if (!("splm" %in% class(x))) stop("'x' must be an 'splm' object.")
      theta(x) <- theta  # Update parameter vector.
   }

   # Calculate mean and standard error of model:
   mu <- predict(x)
   sigma <- predict(x, type = "sigma")

   if (!discrete){
      # Calculate log-probability densities:
      v <- dnorm(x$y, mu, sigma, log = TRUE)
   }else{
      p <- precision(x$y) # Precision of observations.
      yd <- round(x$y/p)
      yd <- yd - min(yd) + 1

      v <- log(pnorm(x$y + p / 2, mu, sigma) - pnorm(x$y - p / 2, mu, sigma))
   }

   return(-sum(v))
}

# THETA<-.SPLM - Parameter assignement method for an 'splm' object.
"theta<-.splm" <- function(x, value, ...) theta(x$theta, ...) <- value

plot.splm <- function(x, xlim, ylim, xlab = "x", ylab = "y", add = FALSE,
                      pch = 21, bg = "grey", col = "red", lwd = 2, data = TRUE, ...){
   # PLOT.SPLM - Plot an 'splm' object.

   # Define plot limits:
   if (missing(xlim)){
      xlim <- NULL;
      if ("x" %in% names(x)) xlim <- range(x$x)
   }
   if (missing(ylim)){
      ylim <- NULL;
      if ("y" %in% names(x)) ylim <- range(x$y)
   }
   if (is.null(xlim)){
      t <- theta(x, "transition", transform = TRUE)
      if (x$order > 1){
         w <- diff(range(t)) / (x$order-1)
         xlim <- c(min(t)-w/2, max(t)+w/2)
      }else{
         w <- theta(x, "window", transform = TRUE)
         xlim <- range(c(t - 2*w, t + 2*w))
      }
   }
   if (is.null(ylim)){
      t <- theta(x, "transition", transform = TRUE)
      t <- sort(unique(c(xlim, t)))
      ylim <- range(predict(x, t))
   }

   # Create an empty plot:
   if (!add) plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, cex.lab = 1.5, ...)

   # Plot data:
   if (data & (all(c("x", "y") %in% names(x)))) points(x$x, x$y, pch = pch, bg = bg, ...)

   # Draw model:
   t <- seq(xlim[1], xlim[2], len = 1000)
   mu <- predict(x, t)
   sigma <- predict(x, t, type = "sigma")
   lines(t, mu, col = col, lwd = lwd, ...)
   lines(t, mu - sigma, lwd = lwd/2, col = col, lty = "dashed")
   lines(t, mu + sigma, lwd = lwd/2, col = col, lty = "dashed")

   # Draw transition points:
   xp <- theta(x, "transition", transform = TRUE)
   #w  <- theta(x, "window", transform = TRUE)
   #if (length(w) == 1) w <- rep(w, length(xp))
   for (i in 1:x$order){
      lines(c(xp[i], xp[i]), c(par("usr")[3], predict(x, xp[i])), lwd = lwd, col = col, lty = "dotted")
   }
}
