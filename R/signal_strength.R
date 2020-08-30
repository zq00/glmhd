#' Estimate problem parameters
#'
#' This function estimates \eqn{(\beta_0, \gamma_0)} using the estimated \eqn{\kappa_s}
#' and intercept \eqn{b}.
#'
#' Assume that \eqn{Y} depends on \eqn{X} as
#' \deqn{
#' \mathrm{P}(Y=1\,|\,X) = \rho'(X^\top \beta + \beta_0),
#' }
#' and let the signal strength be \eqn{\gamma = \mathrm{Var}(X^\top \beta)^{1/2}}.
#' The pair \eqn{(\beta_0, \gamma_0)} satisfies
#' \itemize{
#' \item They are on the phase transition curve \eqn{\kappa(\beta_0, \gamma_0)}.
#'     We estimate the dimension \eqn{\kappa(\beta_0, \gamma_0)} in \code{\link{probe_frontier}}.
#' \item The set of parameters \eqn{(\kappa, \beta_0, \gamma_0)} determines the asymptotic
#'     an estimated intercept \eqn{b}, which should be close to the MLE \code{b_hat}.
#' }
#' We use bisection to search for \eqn{\beta_0}, for each \eqn{\beta_0} we compute \eqn{\gamma_0}
#' on the phase transition curve, and compare the corresponding \eqn{b} with the
#' MLE \code{b_hat}. The algorithm terminates when search window is smaller than \code{eps}.
#'
#' @include prox_op.R h_eq.R find_param.R
#' @param rho_prime A function that computes the success probability \eqn{\rho'(t) = \mathrm{P}(Y=1 | X^\top \beta = t)},
#'     here \eqn{\beta} is the coefficient. The default is logistic model.
#' @param f_prime1 A function. Derivative of the loss function when \eqn{Y = 1}.
#'     The default is the derivative of the negative log-likelihood of logistic
#'     regression when \eqn{Y = 1}.
#' @param f_prime0 A function. Derivative of the loss function when \eqn{Y = -1}.
#'     The default is the derivative of the negative log-likelihood of logistic
#'     regression when \eqn{Y = -1}.
#' @param kappa_hat Numeric. Estimated dimension where the data becomes linearly separable
#' @param kappa Numeric. Problem dimension \eqn{\kappa = p/n}.
#' @param intercept \code{TRUE} if the model contains an intercept.
#' @param b_hat Numeric. MLE of the intercept.
#' @param verbose Should progress be printed? If \code{TRUE}, prints progress at each step.
#' @param eps Numeric. Terminate the algorithm if the search interval is smaller than \code{eps}.
#' @return A list of two components
#' \describe{
#'   \item{gamma_hat}{Estimated signal strength}
#'   \item{beta_hat}{Estimated intercept. \code{NULL} if the model does not have an intercept term}
#' }
#' @examples
#' \dontrun{
#' # no signal case
#' # should return 0, returns 0.0127
#' signal <- signal_strength(kappa_hat = 0.5, intercept = FALSE)
#' signal$gamma_hat
#' signal$b_hat
#' }
signal_strength <- function(rho_prime = rho_prime_logistic,
                            f_prime1 = f_prime1_logistic,
                            f_prime0 = f_prime0_logistic,
                            kappa_hat,
                            kappa = NULL,
                            intercept,
                            b_hat = NULL,
                            verbose = FALSE,
                            eps = 0.1){
  if(verbose) cat("------ Finding signal strength and the intercept ------ \n ")
  if(!intercept){ # no intercept in model
    gamma_hat <- solve_gamma(rho_prime, kappa_hat, 0, FALSE)

    if(verbose) {
      cat("No intercept, estimated signal strength gamma_hat = ", gamma_hat, "\n")
      cat("------ \n")
    }
    return(list(gamma_hat = gamma_hat,
           b_hat = NULL))
  }
  # Search interval
  l <- 0
  # pick the upper bound
  u <- solve_beta(rho_prime, kappa_hat, gamma0 = 0, verbose = FALSE)
  if(u > beta_hat){
    error("Fitted intercept is larger than threshold.")
  }
  if(verbose) cat("Search intercept beta0 in an interval:\n")
  f <- function(beta_hat){
    gamma_hat <- solve_gamma(rho_prime, kappa_hat, beta_hat)
    if(gamma_hat < 0.001) gamma_hat <- 0
    param_new <- find_param(rho_prime, f_prime1, f_prime0,
                            kappa = kappa, gamma = gamma_hat,
                            beta0 = beta_hat, verbose = FALSE)
    # val <- sqrt((sqrt(param_new[1]^2 * gamma_hat^2 + param_new[3]^2) - a_hat)^2 + (param_new[4] - b_hat)^2)
    ifelse(gamma_hat == 0, param_new[3] - b_hat, param_new[4] - b_hat)
  }
  beta_hat <- binary_solve(f, interval = c(l, u), eps = eps, verbose = verbose)
  gamma_hat <- solve_gamma(rho_prime, kappa_hat, beta_hat)
  if(verbose){
    cat("Estimated signal strength gamma_hat = ", gamma_hat, ", intercept beta_0 = ", beta_hat, "\n")
    cat("------ \n")
  }

  list(
    gamma_hat = gamma_hat,
    beta_hat = beta_hat
  )
}


