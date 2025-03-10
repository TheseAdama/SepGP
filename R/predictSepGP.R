#' predictSepGP: Make Predictions Using a SepGP Model
#'
#' @description Compute predictions for new spatial points using a Gaussian Process model with a separable covariance structure (`SepGP`).
#' It computes the posterior mean and covariance matrix for the specified design points based on the model’s prior and covariance functions.
#'
#' @param object An instance of class `SepGP`, representing the Gaussian Process model fitted to training data.
#' @param Dpred Matrix containing new design points for which predictions are to be made. Each row represents a new input point in the same dimension as `D` in `SepGP`.
#' @param covcompute A logical value indicating whether the posterior covariance should be computed. The default value is `TRUE`.
#' @return A list with the following elements:
#'   \item{Mean}{Posterior mean vector of the predictions for each new input in `Dpred`.}
#'   \item{Covsp}{Spatial posterior covariance matrix of the predictions.}
#'   \item{Covmat}{Posterior covariance matrix of the predictions, quantifying the uncertainty for each point in `Dpred`.}
#'
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
#' gp <- SepGP(D, FD, s = 2, covtype = "exp", nugget=1e-12)
#' # Calculate predictions for new design points
#' Dpred <- matrix(seq(0.1, 0.9, length = 5), ncol = 1)
#' Yp <- predictSepGP(gp, Dpred)
#' print(Yp$Mean)
#' print(Yp$Covmat)
setGeneric("predictSepGP", function(object, Dpred, covcompute=TRUE) standardGeneric("predictSepGP"))

#' @rdname predictSepGP
setMethod("predictSepGP", "SepGP", function(object, Dpred, covcompute) {
  V = object@FD - PriorMean(object, object@D)$M

  Sigobs <- Covmat(object, object@D, object@D) + diag(object@nugget, nrow(object@D))

  invSigobs = solve(Sigobs)
  KDDpred <- Covmat(object, object@D, Dpred)
  Mean = PriorMean(object, Dpred)$M + t(KDDpred)%*%invSigobs%*%V

  if(covcompute){
    KDprior = Covmat(object, Dpred, Dpred)
    KDpred = KDprior - t(KDDpred)%*%invSigobs%*%KDDpred
    CC = kronecker(KDpred, object@Sigmat)
  }else{
    CC = KDpred = NULL
  }


  return(list(Mean=Mean, Covsp = KDpred, Covmat=CC))
})
