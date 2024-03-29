% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splm.R
\name{splm}
\alias{splm}
\alias{splm.default}
\title{Smoothed Piecewise Linear Models}
\usage{
splm(x, ...)

\method{splm}{default}(x, theta, weights = FALSE, ...)
}
\arguments{
\item{x}{Numeric vector at which the function is to be evaluated.}

\item{theta}{Named parameter vector with an intercept parameter \code{alpha}, slope parameters \code{beta0, ..., betak},
where \code{k} is the number of transition (i.e. break) points, \code{transition1, ..., transitionk} and
\code{window}, the transition window scale parameter (log-scale). \code{window} may also vary by transition
point.}
}
\value{
If \code{x} is specified, a vector of function evaluations is returned, other wise the function itself is returned.
}
\description{
Functions for defining and evaluating smoothed piecewise-linear model.
}
\section{Methods (by class)}{
\itemize{
\item \code{splm(default)}: Define or evaluate a smoothed piecewise-linear model.

}}
\examples{
# Single transition:
theta <- c(alpha = 0, slope = c(0, 1), transition = 0, window = -3)
x <- seq(-1, 1, length = 1000)
plot(x, splm(theta = theta)(x), type = "l", lwd = 2)
plot(x, splm(theta = theta, weights = TRUE)(x), type = "l", lwd = 2)

# Multiple transitions:
theta <- c(alpha = 0, slope = c(0, 1, -1, 0), transition = c(0, 1, 2), window = -3)
x <- seq(-1, 3, length = 1000)
plot(x, splm(theta = theta)(x), type = "l", lwd = 2)
}
