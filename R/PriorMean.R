#' @title PriorMean: Calculate the Prior Mean for SepGP
#'
#' @description Computes the prior mean for a new set of design points (`Dnew`) based on an existing `SepGP` model.
#' This function constructs a polynomial prior mean function using the coefficients and degree specified in the `SepGP` object.
#'
#' @param object An object of class `SepGP`, representing a Gaussian process model with separable covariance structure.
#' @param Dnew A matrix representing a new design of experiments (DOE) for which the prior mean will be computed.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{M}}{A matrix representing the prior mean values at each point in `Dnew`.}
#'   \item{\code{H}}{A matrix where each row corresponds to the values of the polynomial basis functions evaluated at a point in `Dnew`.}
#' }
#' @export
#' @examples
#' # Test function
#' fxt <- function(x, t) {sin(2 * pi * t * (x - t))}
#' # Design matrix D
#' D <- matrix(seq(0, 1, length = 6), ncol = 1)
#' # Define evaluation points
#' tt <- seq(0, 1, length = 5)
#' # Compute FD
#' FD <- apply(array(tt), 1, function(t){apply(D, 1, function(x){fxt(x,t)})})
#' # Create a SepGP object
#' gp <- SepGP(D, FD, s = 2, covtype = "gauss")
#' # Calculate prior mean for new design points Dnew
#' PriorMean(gp, D)
setGeneric("PriorMean", function(object, Dnew) standardGeneric("PriorMean"))

#' @rdname PriorMean
setMethod("PriorMean", "SepGP", function(object, Dnew) {
  Mprior <- function(D, B, s) {
    H <- t(apply(D, 1, function(u) {
      hvec <- sapply(array(0:s), function(t) sum(apply(array(u), 1, function(u) u^t)))
      hvec <- matrix(hvec, ncol = 1, nrow = s + 1)
      hvec
    }))
    M <- H %*% B
    return(list(M = M, H = H))
  }
  mat <- Mprior(D = Dnew, B = object@B, s = object@s)
  return(mat)
})
