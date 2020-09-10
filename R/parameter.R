#' Model Parameter Class
#'
#' @name parameters
#'
#' @description Class to build, access and manipulate model parameters.
#'
#' @section
#' \describe{
#'    \item{\code{parameter}}{Generic parameter method.}
#'    \item{\code{fix}}{Generic fixed parameter extraction.}
#'    \item{\code{fix<-}}{Generic fixed parameter assignment method.}
#'    \item{\code{theta}}{Generic free parameter extraction.}
#'    \item{\code{theta<-}}{Generic free parameter assignment method.}
#' }
#'

#' @examples
#' p <- parameter(n = 5)
#' p <- parameter(n = 5, label = c("beta0", "beta1", "sigma0", "sigma1", "sigma2"))
#' p <- parameter(n = 5)
#' names(p) <- c("beta0", "beta1", "sigma0", "sigma1", "sigma2")
#' p <- parameter(1:5, label = "sigma")
#' p <- parameter(1:5, label = c("beta0", "beta1", "sigma0", "sigma1", "sigma2"))
#'
#' link(p) <- log
#' link(p, 1:3) <- log
#'
#' # Fix various parameters:
#' fix(p) <- "sigma"
#' fix(p) <- FALSE
#' fix(p) <- TRUE
#' fix(p) <- NULL
#' fix(p) <- c(1, 3, 5)
#' fix(p) <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
#'
#' p <- parameter(1:5, label = c("beta0", "beta1", "sigma0", "sigma1", "sigma2"))
#' fix(p) <- "sigma"
#' link(p, "beta") <- "logit"
#' theta(p) <- c(2.5, 1.5)
#' theta(p, transform = TRUE) <- c(0.25, 0.75)
#'
#' theta(p)
#' theta(p, transform = TRUE)
#' theta(p, full = TRUE)
#' theta(p, transform = TRUE, full = TRUE)
#' theta(p, "sigma", transform = TRUE)
#' fix(p) <- FALSE
#' theta(p, "sigma", transform = TRUE)
#'
#' x <- parameter(1:2, label = "beta")
#' fix(x) <- "beta0"
#' y <- parameter(1:3, label = "sigma")
#' fix(y) <- "sigma"
#' link(y) <- log
#' z <- c(x, y)
#'
#' # 'splm'-type parameter vector:
#' z <- c(parameter(-2, label = "intercept"),
#'        parameter(c(1, -1, 0, 1), label = "slope"),
#'        parameter(n = 3, label = "transition", link = "ordered"),
#'        parameter(0, label = "sigma", link = "log"))
#' theta(z, "transition", transform = TRUE) <- c(0, 1, 2)
#'
#' print(z, format = "short")
#' print(z, format = "long")

c.default <- function(...){
   # C.DEFAULT - Deafult concatenation function.

   # Extract function arguments:
   args <- list(...)

   # Apply usual 'c' function:
   if (!("parameter" %in% unlist(lapply(args, class)))) return(.Primitive("c")(...))

   # Concatenate 'parameter' objects:
   v <- args[[1]]
   if (length(args) > 1) for (i in  2:length(args)) v <- c.parameter(v, args[[i]])

   return(v)
}

c.parameter <- function(x, ...){
   # C.DEFAULT - Deafult concatenation function.

   # Extract function arguments:
   args <- list(...)


   # Extract 'parameter' objects:
   args <- args[unlist(lapply(args, function(x) return("parameter" %in% class(x))))]

   # Single 'parameter' object:
   if (length(args) == 0) return(x)

   # Concatenate multiple objects:
   if (length(args) > 1){
      for (i in 1:length(args)){
         x <- c(x, args[[i]])
      }
      return(x)
   }

   # Extract single parameter object:
   y <- args[[1]]

   if (("parameter" %in% class(x)) & ("parameter" %in% class(y))){
      # Concatenate parameters:
      v <- c(as.numeric.parameter(x), as.numeric.parameter(y))

      # Concatenate fixed parameter vector:
      fixed <- NULL
      if (is.null(attr(x, "fixed"))) fixed <- c(fixed, rep(FALSE, length(x))) else fixed <- c(fixed, attr(x, "fixed"))
      if (is.null(attr(y, "fixed"))) fixed <- c(fixed, rep(FALSE, length(y))) else fixed <- c(fixed, attr(y, "fixed"))
      if (!is.null(fixed)) names(fixed) <- names(v)

      # Concatenate link transforms:
      lx <- link(x)
      ly <- link(y)
      if (!is.null(lx)){
         for (i in 1:length(lx)) lx[[i]]$index <- c(lx[[i]]$index, rep(FALSE, length(v)-length(x)))
         names(lx[[i]]$index) <- names(v)
      }
      if (!is.null(ly)){
         for (i in 1:length(ly)) ly[[i]]$index <- c(rep(FALSE, length(v)-length(y)), ly[[i]]$index)
         names(ly[[i]]$index) <- names(v)
      }
      lv <- c(lx, ly)

      # Concatenate uinverse link transforms:
      lx <- invlink(x)
      ly <- invlink(y)
      if (!is.null(lx)){
         for (i in 1:length(lx)) lx[[i]]$index <- c(lx[[i]]$index, rep(FALSE, length(v)-length(x)))
         names(lx[[i]]$index) <- names(v)
      }
      if (!is.null(ly)){
         for (i in 1:length(ly)) ly[[i]]$index <- c(rep(FALSE, length(v)-length(y)), ly[[i]]$index)
         names(ly[[i]]$index) <- names(v)
      }
      ilv <- c(lx, ly)

      # Concatenate inverse link transforms:
      v <- parameter(v)
      attributes(v) <- c(attributes(v), list(fixed = fixed, link = lv, invlink = ilv))
   }

   return(v)
}


# PARAMETER - Default 'parameter' method.
# 'x' : 'parameter' object or numeric vector.
# 'n' : Number of parameters to be created.
# 'label' : Character string to be used as parameter label.
parameter.default <- function(x, n, label, start = 0, link, drop0 = TRUE, ...){
   # Parse 'n' argument:
   if (!missing(n)){
      if ((length(n) != 1) | !is.numeric(n)) stop("'n' must be a psoitive integer.")
      if ((n < 0) | (n %% 1 != 0)) stop("'n' must be a non-negative integer.")
      if (!missing(x))
         if (length(x) != n) stop("Length of 'x' and specified 'n' are inconsistent.")
   }

   # Parse 'x' argument:
   if (!missing(x)){
      if (!is.numeric(x) | !is.null(dim(x))) stop("'x' must be a numeric vector.")
      if (missing(n)) n <- length(x)
   }

   # Default parameter vector:
   if (missing(x) & missing(n)) x <- as.numeric(NULL)

   # Define empty vector:
   if (missing(x) & !missing(n)) x <- as.numeric(rep(NA, n))

   # Parse 'label' argument:
   if (!missing(label)){
      label <- as.character(label)
      if ((length(label) == 1) & (length(x) > 1) & drop0) label <- paste0(label, start:(n+start-1))
      if (length(label) != length(x))
         stop("'label' length is inconsistent with parameter vector length.")
      names(x) <- label
   }

   # Define class:
   class(x) <- unique("parameter", class(x))

   # Define link function:
   if (!missing(link)) link(x) <- link

   return(x)
}

# SUBSET.PARAMETER - Return a subset of a parameter vector.
subset.parameter <- function(x, ...) return(x[index.parameter(x, ...)])

# INVLINK.PARAMETER - Return the inverse link function(s) associated with a 'parameter' object.
invlink.parameter <- function(x, ...) return(attr(x, "invlink"))

# AS.NUMERIC.PARAMETER - Convert 'parameter object to pure numeric.
as.numeric.parameter <- function(x){
   class(x) <- "numeric"

   tmp <- attributes(x)
   tmp <- tmp[names(tmp) %in% c("names")]
   attributes(x) <- tmp

   return(x)
}

# LINK.PARAMETER - Return the link function(s) associated with a 'parameter' object.
link.parameter <- function(x, ...) return(attr(x, "link"))

# INDEX.PARAMETER - Convert index to a logical vector.
index.parameter <- function(x, index, ...){
   # Reset parameter vector:
   if (missing(index)){
      index <- rep(TRUE, length(x))
      names(index) <- names(x)
   }

   # Index is a character vector of parameter names:
   if (is.character(index)){
      if (all(index %in% names(x))){
         index <- match(index, names(x))
      }else{
         tmp <- index
         index <- grep(index[1], names(x))
         if (length(tmp) > 1){
            for (i in 1:length(tmp)){
               index <- c(index, grep(tmp[i], names(x)))
            }
            index <- sort(unique(index))
         }
      }
      if (length(index) == 0) stop("Unable to match parameter names.")
   }

   # Convert index from numeric to logical:
   if (is.numeric(index)){
      if (any(index %% 1 != 0) | any(index < 1) | any(index > length(x)))
         stop("Invalid index vector for the parameter vector.")
      tmp <- rep(FALSE, length(x))
      tmp[index] <- TRUE
      index <- tmp
   }

   # Check that value is logical:
   if (!is.logical(index)){
      stop("Invalid index for link function specification.")
   }else{
      if (length(index) != length(x))
         stop("Logical index for link function and parameter vector are inconsistent.")
   }

   return(index)
}

# FIX<- Fix target values 'fix' 'parameter' assignment method.
# Usage:
#   fix(x) <- NULL
#   fix(x) <- TRUE
#   fix(x) <- FALSE
#   fix(x) <- index
#   fix(x) <- pname
"fix<-.parameter" <- function(x, value){
   if (is.null(value)){
      # Free all parameters:
      v <- rep(FALSE, length(x))
      names(v) <- names(x)
      attr(x, "fixed") <- v
   }else{
      # Initialize parameter vector:
      if (is.null(attr(x, "fixed"))) fix(x) <- NULL

      # Expand single logical value to whole parameter vector:
      if ((length(value) == 1) & is.logical(value)){
         if (!value) fix(x) <- NULL else attr(x, "fixed") <- (attr(x, "fixed") | TRUE)
         return(x)
      }

      # Extract index of target parameters:
      index <- index.parameter(x, index = value)

      v <- attr(x, "fixed")
      v[index] <- TRUE
      attr(x, "fixed") <- v
   }

   return(x)
}

# LINK<-.PARAMETER - Link function assignement method for an 'parameter' object.
"link<-.parameter" <- function(x, value, ..., add = FALSE){
   # Get logical index:
   index <- index.parameter(x, ...)

   # Extract function string name:
   if (!is.character(value)) value.str <- as.character(substitute(value)) else value.str <- value

   # Determine link inverses for special function values:
   inv.value <- NULL

   # Logarithmic link:
   if (value.str == "log"){
      value <- log
      inv.value <- exp
   }

   # Logit link:
   if (value.str == "logit"){
      inv.value <- function(x) return(exp(x) / (1 + exp(x)))
      value <- function(p) return(log(p) - log(1-p))
   }

   # Ordered vector link:
   if (value.str == "ordered"){
      value <- function(x){
         if (length(x) == 1){
            return(x)
         }else{
            v <- c(x[1], log(diff(x)))
            names(v) <- names(x)
         }
         return(v)
      }
      # Convert to transformed:
      inv.value <- function(x){
         if (length(x) == 1){
            return(x)
         }else{
            v <- cumsum(c(x[1], exp(x[2:length(x)])))
            names(v) <- names(x)
         }
         return(v)
      }
   }

   # Proportion vector link:
   if (value.str == "proportion"){
      value <- function(x){
         p <- 1- sum(x)
         v <- log(x / p)
         names(v) <- names(x)
         return(v)
      }
      # Convert to transformed:
      inv.value <- function(x){
         v <- exp(x) / (1 + sum(exp(x)))
         names(v) <- names(x)
         return(v)
      }
   }

   # Check that 'value' is a function:
   if (!("function" %in% class(value))) stop("'value' must be a function.")

   # Assign link function and index as attributes:
   if (!add | is.null(attr(x, "link"))){
      attr(x, "link") <- list(list(fun = value, index = index, label = value.str))
      if (!is.null(inv.value)) attr(x, "invlink") <- list(list(fun = inv.value, index = index))
   }else{
      attr(x, "link") <- c(attr(x, "link"), list(list(fun = value, index = index, label = value.str)))
      attr(x, "invlink") <- c(attr(x, "invlink"), list(list(fun = inv.value, index = index)))
   }

   return(x)
}

#' @export parameter # PARAMETER - Generic 'parameter' method.
parameter <- function(x, ...) UseMethod("parameter")

#' @export fix
fix        <- function(x, ...) UseMethod("fix")

#' @export "fix<-"
"fix<-"    <- function(x, ...) UseMethod("fix<-")

#' @export theta
theta      <- function(x, ...) UseMethod("theta")

#' @export "theta<-"
"theta<-"  <- function(x, ...) UseMethod("theta<-")

#============================================================================
# THETA.PARAMETER - Extract parameter vector or subsets thereof.
#============================================================================
# Usage:
#    theta(x)        # Returns the set of free parameters.
#    theta(x, index) # Returns a subset of parameters.
#    theta(x, transform = TRUE) # Returns the full transformed parameter vector.
#    theta(x, index, transform = TRUE) # Returns a subset transformed parameter vector. Note that the index applies
#                                      # to the transformed vector and not to the original vector.
#============================================================================
theta.parameter <- function(x, ..., transform = FALSE, full = FALSE){
   # Extract full parameter vector:
   v <- as.numeric.parameter(x)

   # Transform parameter vector:
   if (transform){
      # Transform parameters:
      t <- invlink(x)
      if (!is.null(t)){
         for (i in 1:length(t)){
            w <- t[[i]]$fun(v[t[[i]]$index])
            if (length(w) == length(v[t[[i]]$index])) names(w) <- names(v)[t[[i]]$index]
            v[t[[i]]$index] <- w
         }
      }
   }

   # Remove fixed parameters:
   index <- index.parameter(x, ...)
   if (full) index <- index.parameter(x)

   # Return free parameters:
   if (!full & !is.null(attr(x, "fixed"))){
      index <- index & !attr(x,"fixed")
   }

   return(v[index])
}

#============================================================================
# THETA<-.PARAMETER - Parameter assignement method for a 'parameter' object.
#============================================================================
# Usage:
#    theta(x) <- value            # Assign values to all free parameters.
#    theta(x, index) <- value     # Assign values to subset of parameters.
#    theta(x, transform) <- value # Assign values to complete set of transformed parameters.
#    theta(x, index, transform) <- value # Assign values to subset of transformed parameters.
#============================================================================
"theta<-.parameter" <- function(x, value, ..., transform = FALSE){
   # Extract parameter vector:
   v <- theta(x, transform = transform, full = TRUE)

   # Get subset of parameters:
   if (length(list(...)) == 0){
      # Checked for named parameters:
      if (!is.null(names(value))){
         if (all(names(value) != "")){
            if (all(names(value) %in% names(x))){
               ii <- match(names(value), names(x))
               value <- value[order(ii)]
               index <- index.parameter(x, ii)
            }
         }
      }
   }else{
      index <- index.parameter(x, ...)
   }



   if (!is.null(attr(x, "fixed"))) index <- index & !attr(x, "fixed")

   # Assign value:
   v[index]<- value

   # Assign tranformed values:
   if (transform){
      # Transform parameters:
      t <- link(x)
      if (!is.null(t)){
         for (i in 1:length(t)){
            w <- t[[i]]$fun(v[t[[i]]$index])
            if (length(w) == length(v[t[[i]]$index])) names(w) <- names(v)[t[[i]]$index]
            v[t[[i]]$index] <- w
         }
      }
   }

   # Assign to 'parameter' object:
   attributes(v) <- attributes(x)

   return(v)
}

# PRINT.PARAMETER - Display 'parameter object.
print.parameter <- function(x, format = "long"){
   # Function to generate properly formatted strings:
   bracket <- function(v, quote = "'", sort = FALSE, digits = 4, index){
      if (sort & is.numeric(v)) if (!any(is.na(v))) v <- sort(v)
      if (is.numeric(v)) v <- signif(v, digits)
      if (length(v) > 10){
         n <- length(v)
         if (is.character(v)) v <- paste0(quote, v[c(1, length(v))], quote)
         v <- paste("[", v[1], ", ..., ", v[length(v)], "]", sep = "")
      }else{
         if (is.character(v)) v <- paste(quote, v, quote, sep = "")
         if (!missing(index) & is.numeric(v)){
            v <- as.character(v)
            v[index] <- paste0("(", v[index], ")")
         }
         if (length(v) > 1) v <- paste0("[", paste(v, collapse = ", "), "]")
      }
      return(v)
   }

   # Display summary of 'alkey' object:
   cat(paste0("'parameter' object:\n"))

   # Determine parameter groupings if long-form is desired:
   flag <- FALSE
   if (format == "long"){
      # Long-form parameter display:
      v <- theta(x, full = TRUE)
      if (!is.null(names(v))){
         if (all(names(v) != "")){
            if (all(gsub("^[a-zA-Z_.]+[0-9]*$", "", names(v)) == "")){
               group <- gsub("[0-9]*$", "", names(v))
               groups <- unique(group)
               if (length(groups) > 1){
                  flag <- TRUE
               }
            }
         }
      }
   }

   if ((format == "short") | !flag) {
      # Short-form parameter display:
      cat(paste0("            theta : ", bracket(theta(x, full = TRUE), index = attr(x, "fixed")), "\n"))
      cat(paste0("      transformed : ", bracket(theta(x, transform = TRUE, full = TRUE), index = attr(x, "fixed")), "\n"))
   }else{
      # Display parameter vector:
      for (i in 1:length(groups)){
         if (i == 1) str <- "            theta " else str <- "                  "
         index <- group == groups[i]
         cat(paste0(str, ": ", bracket(groups[i]), " = ",  bracket(v[index], index = attr(x, "fixed")[index]), "\n"))
      }
      v <- theta(x, transform = TRUE, full = TRUE)
      # Display transformed parameter vector:
      for (i in 1:length(groups)){
         if (i == 1) str <- "      transformed " else str <- "                  "
         index <- group == groups[i]
         cat(paste0(str, ": ", bracket(groups[i]), " = ",  bracket(v[index], index = attr(x, "fixed")[index]), "\n"))
      }
   }

   # Display link functions:
   if (!is.null(link(x))){
      t <- link(x)
      cat(paste0("          link(s) : ", bracket(unlist(lapply(t, function(x) x$label))), "\n"))
   }
}
