#' 2-dimensional Gaussian Expectation
#'
#' Expected value of a function of two independent
#' standard Gaussian variables
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
  integrand <- function(x) f(x) * dnorm(x[1]) * dnorm(x[2])

  hcubature(integrand, lowerLimit = c(-8,-8), upperLimit = c(8,8),
                      fDim = 1, maxEval = 0, tol = 1e-4,  ...)$integral
}
