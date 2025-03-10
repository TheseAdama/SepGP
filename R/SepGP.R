#' @title SepGP: Gaussian Processes with Separable Covariance Functions
#'
#' @description The `SepGP` class implements Gaussian processes with separable covariance functions for regression and interpolation tasks.
#' It provides flexibility in specifying covariance structures, priors, and additional regularization through the nugget effect, allowing efficient hyperparameter optimization.
#'
#' @slot D Design of experiments matrix containing the input points.
#' @slot FD Matrix representing evaluations of the computer code at the points in `D`.
#' @slot s Degree of the polynomial used for the prior mean function.
#' @slot B Matrix of coefficients for the prior mean, based on `D`.
#' @slot covtype Type of covariance function; options include `"gauss"`, `"exp"`, `"matern3_2"`, and `"matern5_2"`.
#' @slot phi Vector of hyperparameters controlling the covariance structure.
#' @slot philow Lower bounds for the hyperparameters in `phi`.
#' @slot phiup Upper bounds for the hyperparameters in `phi`.
#' @slot Sigmat Noise covariance matrix, representing observation noise in the data.
#' @slot nugget A small positive value added to the diagonal of the covariance matrix to improve numerical stability, particularly in cases with near-zero variance or when data points are very close to each other. Defaults to `1e-8`.
#' @exportClass SepGP
setClass("SepGP",
         slots = list(D = "matrix",
                      FD = "matrix",
                      s = "numeric",
                      B = "matrix",
                      covtype= "character",
                      phi = "numeric",
                      philow="numeric",
                      phiup = "numeric",
                      Sigmat = "matrix",
                      nugget = "numeric"))

#' Constructor for SepGP Class
#'
#' @description Constructs an instance of the `SepGP` class, initializing a Gaussian process model with customizable covariance, prior mean structures, and optional nugget effect.
#'
#' @param D Matrix representing the design of experiments (input points).
#' @param FD Matrix representing the computer code evaluations on `D`.
#' @param s Degree of the polynomial for the prior mean function.
#' @param B Optional. Coefficients for the prior mean on `D`. If not specified, defaults to a matrix of ones.
#' @param covtype Type of covariance function to be used, with default value `"gauss"`. Acceptable values are `"gauss"`, `"exp"`, `"matern3_2"`, and `"matern5_2"`.
#' @param phi Optional. Hyperparameter vector for the covariance function. If not specified, defaults to a vector of ones.
#' @param philow Optional. Lower bounds for each hyperparameter in `phi`. If not specified, defaults to `0.1` for each element.
#' @param phiup Optional. Upper bounds for each hyperparameter in `phi`. If not specified, defaults to `10` for each element.
#' @param Sigmat Optional. Temporal covariance matrix. If not specified, defaults to an identity matrix of appropriate size.
#' @param nugget Optional. A small positive value added to the diagonal of the covariance matrix to improve numerical stability, defaulting to `1e-8`.
#'
#' @return An object of class `SepGP` representing a Gaussian process model with the specified settings.
#' @importFrom methods new
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
SepGP <- function(D, FD, s, B = NULL, covtype="gauss", phi = NULL, philow = NULL, phiup = NULL,
                  Sigmat = NULL, nugget=1e-16) {

  # Vérification des arguments
  if (!is.matrix(D)) stop("The design matrix D must be a matrix.")
  if (!is.matrix(FD)) stop("The observations matrix FD must be a matrix.")
  if (nrow(FD) != nrow(D)) stop("The number of rows in FD must be equal to the number of rows in D.")
  if (!is.numeric(s) || length(s) != 1) stop("s must be a single numeric value.")
  if (!is.character(covtype) || !(covtype %in% c("gauss", "exp", "matern3_2", "matern5_2"))) {
    stop("covtype must be one of 'gauss', 'exp', 'matern3_2', or 'matern5_2'.")
  }

  nD <- ncol(D)

  # Valeurs par défaut pour phi, philow, phiup
  if (is.null(phi)) {
    phi <- rep(1, nD)
  } else {
    if (!is.numeric(phi) || length(phi) != nD) stop("phi must be numeric with length equal to the number of columns in D.")
  }

  if (is.null(philow)) {
    philow <- rep(0.1, nD)  # Valeur par défaut pour philow
  } else {
    if (!is.numeric(philow) || length(philow) != nD) stop("philow must be numeric with length equal to the number of columns in D.")
  }

  if (is.null(phiup)) {
    phiup <- rep(10, nD)
  } else {
    if (!is.numeric(phiup) || length(phiup) != nD) stop("phiup must be numeric with length equal to the number of columns in D.")
  }

  if (is.null(B)) {
    # Valeur par défaut pour B
    B <- matrix(1, ncol = ncol(FD), nrow = s + 1)
  } else {
    if (!is.matrix(B)) stop("B must be a matrix.")
    if (ncol(B) != (s + 1)) stop("The number of columns in B must be equal to s + 1.")
    if (nrow(B) != nrow(FD)) stop("The number of rows in B must be equal to the number of rows in FD.")
  }

  if (is.null(Sigmat)) {
    # Valeur par défaut pour Sigmat
    Sigmat <- diag(1, nrow = ncol(FD))
  } else {
    if (!is.matrix(Sigmat)) stop("Sigmat must be a matrix.")
    if (ncol(Sigmat) != ncol(FD)) stop("The size of Sigmat must be equal to the number of rows in FD.")
  }

  if (!is.numeric(nugget) || length(nugget) != 1 || nugget <= 0 || nugget >= 1) {
    stop("The nugget must be a positive numeric scalar less than 1.")
  }

  # Création de l'objet SepGP
  new("SepGP",
      D = D,
      FD = FD,
      s = s,
      B = B,
      covtype = covtype,
      phi = phi,
      philow = philow,
      phiup = phiup,
      Sigmat = Sigmat,
      nugget=nugget)
}



