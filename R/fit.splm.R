fit.splm <- function(x, trace = FALSE, print = TRUE, niter = 3, ...){
   # FIT.SPLM - Fit an 'splm' object to data.
   
   # Initialize parameter vector:
   if (any(is.na(theta(x, full = TRUE)))) x <- init(x)
   
   # Store starting log-likelihood value:
   ll <- log.likelihood.splm(x, ...)
   
   if (!trace) trace <- 0
   for (i in 1:niter){
      v <- optim(theta(x), log.likelihood.splm, control = list(trace = 0, maxit = 5000), x = x, ...)
      theta(x) <- v$par
   }
   
   if (print){
      str <- ""
      if (v$convergence == 0) str <- paste0(str, "Converged") else str <- paste0(str, "Did not converge.") 
      str <- paste0(str, ", ln(theta) = ", as.character(v$value))
      str <- paste0(str, ", ln(delta) = ", as.character(log(ll - v$value)))
      str <- paste0(str, "\n")
      cat(str)
   }
   
   return(x)
}
