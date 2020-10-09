variogram.variogram <- function(x, h){
   # VARIOGRAM - Evaluate variogram at specified distances.
   
   v <- x$vfun(h, nugget = x$nugget, sill = x$sill, range = x$range)
   
   return(v)
}
