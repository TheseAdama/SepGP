#' Covmat: Calculate the Covariance Matrix for SepGP
#'
#' This function computes the covariance matrix between two DOEs based on a SepGP model.
#'
#' @param object An object of class SepGP
#' @param D1 First DOE
#' @param D2 Second DOE
#' @return Covariance matrix
#' @export
#' @examples
#' # Test function
#' fxt <- function(x, t) {sin(2 * pi * t * (x - t))}
#' # Design matrix D
#' D <- matrix(seq(0, 1, length = 10), ncol = 1)
#' # Define evaluation points
#' tt <- seq(0, 1, length = 10)
#' # Compute FD
#' FD <- apply(array(tt), 1, function(t){apply(D, 1, function(x){fxt(x,t)})})
#' # Create a SepGP object
#' gp <- SepGP(D, FD, s = 2, covtype = "gauss")
#' Covmat(gp, D, D)
setGeneric("Covmat", function(object, D1, D2) standardGeneric("Covmat"))

#' @rdname Covmat
setMethod("Covmat", "SepGP", function(object, D1, D2) {
  Covsep <- function(theta1, theta2, phi, covtype) {
    r = ((theta1 - theta2)/phi)
    r = norm(r,"2")

    if(covtype=="exp"){o = exp(-r)}

    if(covtype=="gauss"){
      o = exp(-0.5*r**2)}

    if(covtype=="matern3_2"){o = (1+ sqrt(3)*r)*exp(-sqrt(3)*r)}

    if(covtype=="matern5_2"){
      o = (1+ sqrt(5)*r + (5/3)*r^2)*exp(-sqrt(5)*r)}

    return(o)
  }

  mat <- outer(1:nrow(D1), 1:nrow(D2), Vectorize(function(i, j) { Covsep(D1[i, ], D2[j, ], phi=object@phi, covtype=object@covtype)}))
  return(mat)
})
