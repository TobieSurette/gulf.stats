#' Generic class functions:
#'
#' @name gulf.stats
#'
#' @description \code{gulf.stats} generic class functions.
#'
#' @param x Object.
#'
#' @section Generic methods:
#' \describe{
#'    \item{\code{init}}{Generic method for initializing objects.}
#'    \item{\code{loglike}}{Generic log-likelihood method.}
#'    \item{\code{fit}}{Generic model fitting method.}
#'    \item{\code{link}}{Generic link function method.}
#'    \item{\code{prior}}{Generic \code{prior} extraction method.}
#'    \item{\code{prior<-}}{Generic \code{prior} assignment method.}
#' }
#'

#' @rdname gulf.stats
#' @export init
init       <- function(x, ...) UseMethod("init")

#' @rdname gulf.stats
#' @export loglike
loglike    <- function(x, ...) UseMethod("loglike")

#' @rdname gulf.stats
#' @export fit
fit        <- function(x, ...) UseMethod("fit")

#' @rdname gulf.stats
#' @export link
link       <- function(x, ...) UseMethod("link")

#' @rdname gulf.stats
#' @export "link<-"
"link<-"   <- function(x, ...) UseMethod("link<-")

#' @rdname gulf.stats
#' @export prior
"prior"    <- function(x, ...) UseMethod("prior")

#' @rdname gulf.stats
#' @export "prior<-"
"prior<-"  <- function(x, ...) UseMethod("prior<-")


