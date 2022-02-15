#' 2-dimensional Gaussian Expectation
#'
#' Expected value of a function of two independent
#' standard Gaussian variables
#'
#' @return A numeric value
#' \deqn{
#' \int \int f(x, y) \varphi(x) \varphi(y) dx dy,
#' }
#' where \eqn{ \varphi(\cdot)} is the standard Gaussian density.
#'
#' @param f Function to integrate over. Input to
#' \code{f} should be a length 2 vector.
#' @param ... Additional arguments to \link[cubature]{hcubature}
#' @importFrom cubature hcubature
#' @importFrom stats dnorm
#' @examples
#' \dontrun{
#' # Second moment of a Gaussian
#' # tol = 0.0001, true value = 2, evaluated = 2.000001
#' f <- function(x) x[1]^2 + x[2]^2
#' integrate2_normal(f, tol = 1e-4)
#' }

integrate2_normal <- function(f, ...){
  integrand <- function(x){
    matrix(apply(x, 2, function(t) f(t) * exp(-sum(t^2) / 2) / 2 / pi), ncol = ncol(x))
  }

  pcubature(integrand, lowerLimit = c(-10,-10), upperLimit = c(10,10),
            fDim = 1, maxEval = 0, absError = 0, vectorInterface = T)$integral
}
