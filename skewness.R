skewness <- function(x){
   # SKEWNESS - Calculate sample skewness:
   
   # Number of observations:
   n <- length(x)
   
   # Degenerate values:
   if (n < 3) return(NA)
   s <- sd(x)
   if (s == 0) return(0)

   # Calculate sample skewness:
   v <- n /((n-1)*(n-2)) * sum(((x-mean(x))/s)^3)

   return(v)
}
