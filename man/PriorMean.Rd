% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PriorMean.R
\name{PriorMean}
\alias{PriorMean}
\alias{PriorMean,SepGP-method}
\title{PriorMean: Calculate the Prior Mean for SepGP}
\usage{
PriorMean(object, Dnew)

\S4method{PriorMean}{SepGP}(object, Dnew)
}
\arguments{
\item{object}{An object of class \code{SepGP}, representing a Gaussian process model with separable covariance structure.}

\item{Dnew}{A matrix representing a new design of experiments (DOE) for which the prior mean will be computed.}
}
\value{
A list with the following elements:
\describe{
\item{\code{M}}{A matrix representing the prior mean values at each point in \code{Dnew}.}
\item{\code{H}}{A matrix where each row corresponds to the values of the polynomial basis functions evaluated at a point in \code{Dnew}.}
}
}
\description{
Computes the prior mean for a new set of design points (\code{Dnew}) based on an existing \code{SepGP} model.
This function constructs a polynomial prior mean function using the coefficients and degree specified in the \code{SepGP} object.
}
\examples{
# Test function
fxt <- function(x, t) {sin(2 * pi * t * (x - t))}
# Design matrix D
D <- matrix(seq(0, 1, length = 6), ncol = 1)
# Define evaluation points
tt <- seq(0, 1, length = 5)
# Compute FD
FD <- apply(array(tt), 1, function(t){apply(D, 1, function(x){fxt(x,t)})})
# Create a SepGP object
gp <- SepGP(D, FD, s = 2, covtype = "gauss")
# Calculate prior mean for new design points Dnew
PriorMean(gp, D)
}
