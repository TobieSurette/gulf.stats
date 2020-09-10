"fix<-.splm" <- function(x, value, ...){
   # FIX<- 'fix' 'splm' assignment method.
   
   fix(x$theta) <- value
      
   return(x)
}
