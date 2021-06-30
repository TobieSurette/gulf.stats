#' @title Polynomial Class
#'
#' @description Methods for creating, evaluating, manupulating and displaying polynomial functions.
#'
#' @param beta Numeric vector of polynomial coefficients. First element is the constant coefficient,
#'             the second element is the linear coefficient, the third is the quadratic coefficient, etc...
#' @param x,y  Polynomial object.
#'
#' @examples
#' p1 <- polynomial(1:3)    # Create polynomial 1 + 2*x + 3*x^2.
#' p2 <- polynomial(0:5)    # Create polynomial x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5.
#'
#' # Evaluation:
#' x <- 1:10
#' polynomial(1)(x)        # Constant function.
#' polynomial(c(1,0,2))(x) # Quadratic 1 + 2*x^2.
#' p1(x)                   # Evaluate 'p1' at 'x'.
#' (p1 + p2)(x)            # Evaluate the sum of 'p1' and 'p2' at 'x'.
#' (p1 - p2)(x)            # Evaluate the difference of 'p1' and 'p2' at 'x'.

#' @export
polynomial <- function(beta) UseMethod("polynomial")

#' @describeIn polynomial Polynomial object creation and evaluation.
#' @export
polynomial.default <- function(beta){
   if ("polynomial" %in% class(beta)) return(beta)

   # Simplifify polynomial:
   if (all(beta == 0)) beta <- 0

   # Horner's method for polynomial evaluation:
   f <- function(x){
      y <- beta[[length(beta)]] * rep(1, length(x))
      if (length(beta) > 1) for (i in 1:(length(beta)-1)) y <- y * x + beta[length(beta)-i]
      return(y)
   }
   class(f) <- c("polynomial", class(f))
   attr(f, "beta") <- beta

   return(f)
}

#' @describeIn polynomial Polynomial unary and binary addition operator.
#' @export
"+.polynomial" <- function(x, y){
   if (missing(y)) return(x)

   # Convert to polynomial:
   x <- polynomial(x)
   y <- polynomial(y)

   # Order of largest polynomial:
   k <- max(length(coef(x)), length(coef(y)))

   # Add coefficients:
   beta <- c(coef(x), rep(0, k - length(coef(x)))) +
      c(coef(y), rep(0, k - length(coef(y))))

   p <- polynomial(beta)

   return(p)
}

#' @describeIn polynomial Polynomial unary and binary substraction operator.
#' @export
"-.polynomial" <- function(x, y){
   if (missing(y)) return(polynomial(-coef(x)))
   y <- polynomial(-coef(y))

   return(x + y)
}

#' @describeIn polynomial Polynomial multiplication operator.
#' @export
"*.polynomial" <- function(x, y){
   # Convert to polynomial:
   x <- polynomial(x)
   y <- polynomial(y)

   # Extract coefficients:
   a <- coef(x)
   b <- coef(y)

   # Calculate coefficients:
   p <- rep(0, length(a) + length(b))
   for (i in 1:length(a)){
      for (j in 1:length(b)){
         p[i+j-1] <- p[i+j-1] + a[i] * b[j]
      }
   }

   return(polynomial(p))
}

#' @describeIn polynomial Polynomial print method.
#' @export
print.polynomial <- function(x){
   s <- as.character(x)

   # Print result:
   cat(paste0("p(x) = ", s, "\n"))
}

#' @rawNamespace S3method(as.character,polynomial)
as.character.polynomial <- function(x){
   # Prepare polynomial coefficients:
   beta <- coef(x)
   ix <- 0:(length(beta)-1)
   ix <- ix[beta != 0]
   beta <- beta[beta != 0]
   if (length(beta) == 0){
      beta <- 0
      ix <- 0
   }

   # Build formula representation:
   s <- paste0(beta, "x^", ix)
   s <- paste0(s, collapse = " + ")
   s <- gsub("x^0 ", " ", s, fixed = TRUE)
   s <- gsub(" 1x", " x", s, fixed = TRUE)
   s <- gsub("^0 [+] ", "", s)
   s <- gsub("0x[\\^][0-9]+", "", s)
   s <- gsub("x[\\^]1 ", "x ", s)
   if (s == " + ") s <- ""
   if (s == "") s <- "0"
   s <- paste0("p(x) = ", s)
   s <- gsub(" 1x ", " x ", s)
   s <- gsub(" +$", "", s)
   s <- gsub("x[\\^]0$", "", s)
   s <- gsub("^p[(]x[)] = ", "", s)
   s <- gsub("+ -", "- ", s, fixed = TRUE)

   return(s)
}


#' @describeIn polynomial Polynomial print method.
#' @export
coef.polynomial <- function(x) return(attr(x, "beta"))

