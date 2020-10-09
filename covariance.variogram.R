covariance.variogram <- function(x, h){
   # COVARIANCE - Evaluate covariance associated with variogram at specified distances.
   
   v <- x$sill - variogram(x, h)
   
   return(v)
}
