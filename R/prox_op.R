#' One-dimensional proximal operator
#'
#' @param f_prime Function. Derivative of the function f.
#' @param lambda Numeric. Inverse penalty.
#' @param x Numeric. Input variable.
#' @return A numeric value
#'     \deqn{\mathrm{prox}_{\lambda f}(x) = \mathrm{argmin}_x f(z) + \frac{1}{2\lambda}(z-x)^2.}
#'     Solution to the equation
#'     \deqn{f'(z) + \frac{1}{\lambda}(z-x) = 0.}
#' @examples
#' f <- function(x) 1/(1+exp(-x))
#' prox_op(f, 1, 1)

prox_op <- function(f_prime, lambda, x){
  f <- function(z)  z + lambda * f_prime(z) - x
  uniroot(f, interval = c(-100,100), extendInt = "yes")$root
}
