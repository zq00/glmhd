#' Logistic loss function
#'
#' Success probability and derivatives of the
#' negative log-likelihood of a logistic model
#'
#' \code{rho_prime_logistic} computes the successs probability of
#' a logistic model \eqn{\rho'(t) = 1/(1+e^{-t})}.
#'
#' \code{f_prime1_logistic} and \code{f_prime0_logistic} are
#' the derivative of the negative log-likelihood when \eqn{Y=1} and \eqn{Y=-1}.
#'
#' @param x Numeric. Input variable.
#' @export

rho_prime_logistic <- function(x) 1 / (1 + exp(-x))

#' @rdname rho_prime_logistic
#' @export
f_prime1_logistic <- function(x) -1 / (1 + exp(x))

#' @rdname rho_prime_logistic
#' @export
f_prime0_logistic <- function(x) 1 / (1 + exp(-x))


