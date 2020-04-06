#' Line intersections
#' 
#' Returns coordinates where pairs of lines intersect.
#' 
#' @param beta0 Two-element vector or two-column matrix of line slopes and intercept values.
#' @param beta1 Two-element vector or two-column matrix of line slopes and intercept values.
#' @param pairwise Logical value specifying whether to calculate intersection points for all possible slope and intercept pair combinations.
#'  
#' @return Returns a list of x and y intersection points.
#' 
#' 
#' # Example:
#' n <- 10
#' m <- 2*runif(n)-1
#' b <- 2*runif(n)-1
#' plot(c(-2, 2), c(-2, 2), type = "n", xlab = "x", ylab = "y")
#' for (i in 1:n) abline(b[i], m[i], col = "grey70")
#' r <- line.intersect(cbind(m, b))
#' points(r$x, r$y, pch = 21, bg = "red")
#' 
#' @export line.intersect

line.intersect <- function(beta0, beta1, pairwise = TRUE){
   
   # Convert vector to matrix:
   if (is.null(dim(beta0))) beta0 <- t(beta0)
   
   # Duplicate single input:
   if (missing(beta1)){
      beta1 <- beta0
      if (is.null(dim(beta1))) beta1 <- t(beta1)
   } 
   
   # Slopes:
   m0 <- beta0[,1]
   m1 <- beta1[,1]
   
   # Intercepts:
   b0 <- beta0[,2]
   b1 <- beta1[,2]
  
   # Expand input values to all possible pair combination:
   if (pairwise){
      n <- c(length(m0), length(m1))
      m0 <- repvec(m0, ncol = n[2])
      b0 <- repvec(b0, ncol = n[2])
      m1 <- repvec(m1, nrow = n[1])
      b1 <- repvec(b1, nrow = n[1])   
      m0 <- as.numeric(m0); m1 <- as.numeric(m1)
      b0 <- as.numeric(b0); b1 <- as.numeric(b1)
   }
   
   # Intersections:
   x <- (b1 - b0) / (m0 - m1)
   y <- m0 * x + b0
   
   # Remove duplicates and undefined values:
   ii <- !is.na(x) & !is.na(y) & !duplicated(cbind(x, y))
   x <- x[ii]; y <- y[ii]
   
   return(data.frame(x = x, y = y))
}
