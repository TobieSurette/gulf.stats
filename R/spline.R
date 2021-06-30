#' @title Spline Functions
#'
#' @description Spline function class
#'
#' @param polynomial Polynomial coefficients, polynomial object or list of polynomial objects.
#' @param knots Numeric vector of spline knots.
#'

#' @export
spline <- function(x) UseMethod("spline")

#' @export spline.default
spline.default <- function(x, polynomial, knots){
   # Default method:
   if (!missing(x)) if (is.numeric(x)) return(stats::spline(x, ...))

   # Parse input arguments:
   if (!is.numeric(knots)) stop("'knots' must be a numeric vector.")
   if ((length(polynomial)-1) != length(knots)) stop("Number of polynomials must be one more than the number of knots.")

   # Parse polynomial argument:
   p <- polynomial
   if (is.numeric(p)) p <- polynomial(p)
   if ("polynomial" %in% class(p)) p <- list(p)

   # Expand knots:
   knots <- c(-Inf, sort(knots), Inf)

   # Evaluation function:
   f <- function(x){
      k <- length(knots)-1

      y <- rep(NA, length(x))
      for (i in 1:k){
         ix <- (x >= knots[i]) & (x < knots[i+1])
         y[ix] <- p[[i]](x[ix])
      }

      return(y)
   }

   attr(f, "polynomial") <- p
   attr(f, "knots") <- knots[-c(1, length(knots))]
   class(f) <- c("spline", class(f))

   return(f)
}

#' @export
print.spline <- function(x){
   p <- attr(x, "polynomial")
   s <- lapply(p, as.character)
   knots <- c(-Inf, attr(x, "knots")  , Inf)
   s <- unlist(s)
   s <- paste0("   p_", 1:length(s), "(x) = ", s)
   space <- NULL
   for (i in 1:length(s)) space[i] <- paste(rep(" ", max(nchar(s)) - nchar(s[i]) + 3), collapse = "")
   s <- paste0(s, ",", space, " x in [",  knots[-length(knots)], ", ", knots[-1], "]")
   cat("'spline' object:\n")
   for (i in 1:length(s)){
     cat(paste(s[i], "\n"))
   }
   cat("\n")
}


