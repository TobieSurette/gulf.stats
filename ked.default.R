ked.default <- function(x, x0, y, z, z0, nugget, sill, range, positive = TRUE){
   # KED - Perform point kriging with external drift.
   #   'x'  : Data coordinates.
   #   'x0' : Prediction coordinates.
   #   'y'  : Data response values.
   #   'z'  : External drift values at data locations.
   #   'z0' : External drift values at prediction locations.

   # Number of data points:
   n <- nrow(x)  
      
   # Basic functions:
   step <- function(x) return((x >= 0) + 1 - 1)
   spherical <- function(h) return((1-step(h-1)) * (1 - (1.5*h - .5*(h^3))))   
      
   # Nugget value contribution to covariance matrix:
   k <- cbind(diag(rep(nugget, n)), 0) 
   
   # Sill contribution to covariance matrix:   
   t <- as.matrix(rbind(x, x0) / range)  
   t <- t %*% t(t)
   h = sqrt(-2*t+diag(t) %*% matrix(1,ncol = n+1) + matrix(1, nrow = n+1) %*% t(diag(t)))
   h = h[1:n,]
   k <- k + sill * spherical(h)

   # Separate covariance matrix components:
   k0 = k[,n+1] # Prediction versus data covariance values.
   k = k[,1:n]  # Data covariance values.

   # ExtetrAdd external drift variable:
   k <- rbind(cbind(k, 1), 1)
   k[n+1, n+1] <- 0
   k <- rbind(cbind(k, c(z, 0)), c(z, 0, 0))
   k0 <- matrix(c(k0, 1), ncol = 1)
   k0 <- rbind(k0, z0)
   
   # Solve cokriging system by Gaussian elimination:
   const <- max(diag(k)) / max(max(abs(k[,(n+1):ncol(k)])))
   k[1:n,1:n] <- k[1:n,1:n] / const
   k0[1:n,] <- k0[1:n,] / const
   
   # Get kriging weights:
   w <- solve(k) %*% k0  
  
   # calculation of cokriging estimates
   x0s = sum(w[1:n,] * y)
    
   # calculation of cokriging variances
   s = (nugget + sill) / const
   t = t(w) %*% k0
   s = const * (s-t);

   # Set negative values to zero:
   if (positive) x0s[x0s < 0] <- 0
    
   return(list(mean = x0s, var = s, weight = w, K = k, K0 = k0))
}
