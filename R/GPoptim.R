#' GPoptim: Hyperparameter Optimization for a SepGP Model
#'
#' The `GPoptim` function optimizes the hyperparameters (e.g., spatial covariance function parameters) of a given `SepGP` model
#' using a leave-one-out cross-validation method. This optimization seeks to minimize the error in predicting points
#' left out of the training set, enhancing the model's predictive accuracy and generalization.
#'
#' @param object An instance of the `SepGP` class representing the Gaussian process model to be optimized.
#'               The `SepGP` object should already contain the initial design matrix (`D`), observed data (`FD`),
#'               polynomial degree (`s`), covariance function type (`covtype`), and initial hyperparameter values.
#' @param maxit (Optional) Integer specifying the maximum number of iterations for the optimization algorithm (refer to the `optim` function in the `stats` package for details). Default is 1000.
#'
#' @return An optimized `SepGP` object with updated hyperparameters:
#' - `phi`: Optimal values of the covariance hyperparameters.
#' - `B`: Coefficients of the prior mean based on the optimized hyperparameters.
#' - `Sigmat`: Temporal covariance matrix adjusted according to the new hyperparameters.
#'
#' @export
#' @importFrom stats optim
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
#' gp <- GPoptim(gp, maxit=10)
#' print(gp@phi)
#' print(gp@B)
#' print(gp@Sigmat)
setGeneric("GPoptim", function(object, maxit=1000) standardGeneric("GPoptim"))

#' @rdname GPoptim
setMethod("GPoptim", "SepGP", function(object, maxit=1000) {

  if(!is.numeric(maxit)) stop("maxit must be numeric !")

  OLOO <- function(phi){
    m = nrow(object@D)
    E = rep(0, m)
    for(i in 1:m){
      Oi = SepGP(D = matrix(object@D[-i, ], ncol=ncol(object@D)),
                 FD = object@FD[-i,],
                 s = object@s,
                 covtype = object@covtype,
                 phi = phi,
                 philow = object@philow,
                 phiup = object@phiup,
                 nugget = object@nugget)

      H <- t(apply(Oi@D, 1, function(u) {
        hvec <- sapply(array(0:Oi@s), function(t) sum(apply(array(u), 1, function(u) u^t)))
        hvec <- matrix(hvec, ncol = 1, nrow = Oi@s + 1)
        hvec
      }))
      Sigobs <- Covmat(Oi, Oi@D, Oi@D) + diag(Oi@nugget, nrow(Oi@D))
      invSigobs <- solve(Sigobs)
      Oi@B <- (t(H)%*%invSigobs%*%H)%*%(t(H)%*%invSigobs%*%Oi@FD)

      V <- Oi@FD - PriorMean(Oi, Oi@D)$M
      Oi@Sigmat <- (1/nrow(Oi@D)) * t(V)%*%invSigobs%*% V

      R <- predictSepGP(Oi, matrix(object@D[i, ], nrow=1), covcompute = FALSE)
      E[i] <- sum((matrix(object@FD[i, ], nrow=1) - R$Mean)^2)
      #sum(log(diag(R$Covmat)+1e-16) + ((matrix(object@FD[i, ], nrow=1) - R$Mean)^2/diag(R$Covmat)))
    }
    o <- sum(E, na.rm = TRUE)
    return(o)
  }


  R <- stats::optim(par = object@phi,
                    fn = OLOO,
                    lower = object@philow,
                    upper = object@phiup,
                    method = "L-BFGS-B",
                    control = list(maxit = maxit, trace = 1, REPORT = 1))
  phistar = R$par

  object@phi = phistar
  H <- t(apply(object@D, 1, function(u) {
    hvec <- sapply(array(0:object@s), function(t) sum(apply(array(u), 1, function(u) u^t)))
    hvec <- matrix(hvec, ncol = 1, nrow = object@s + 1)
    hvec
  }))
  Sigobs <- Covmat(object, object@D, object@D) + diag(object@nugget, nrow(object@D))
  invSigobs <- solve(Sigobs)
  object@B <- (t(H)%*%invSigobs%*%H)%*%(t(H)%*%invSigobs%*%object@FD)

  V <- object@FD - PriorMean(object, object@D)$M
  object@Sigmat <- (1/nrow(object@D))*t(V)%*%invSigobs%*% V

  return(object)
})

