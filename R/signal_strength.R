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
                            kappa_hat,
                            intercept = FALSE,
                            p0 = NULL,
                            verbose = FALSE,
                            tol = 1e-4){
  if(verbose) cat("-- Finding signal strength and the intercept -- \n ")
  if(!intercept){ # no intercept in model
    gamma_hat <- solve_gamma(rho_prime, kappa_hat, beta0 = 0, verbose = FALSE)
    if(verbose) { cat("No intercept, estimated signal strength gamma_hat = ", gamma_hat, "\n ---- \n")}
    return(list(gamma_hat = gamma_hat))
  }
  # expected proportion of Y= 1
  prop <- function(beta0, gamma){
    f <- function(z){dnorm(z) * rho_prime(beta0 + gamma * z)}
    integrate(f, -8, 8, abs.tol = 1e-5)$value
  }
  # a system of two equations
  f <- function(x){
    beta0 <- x[1]; gamma <- x[2]
    val <- c(kappa_hat - solve_kappa(rho_prime = rho_prime, abs(beta0), gamma),
             p0 - prop(beta0, gamma))
  }
  # initial values
  x_init <- c(log(p0 / (1-p0)), 2)
  sol <- fsolve(f, x_init, J = NULL, maxiter = 100, tol = tol)$x

  if(verbose){cat("Estimated signal strength gamma_hat = ", sol[2], ", intercept beta_0 = ", sol[1], "\n ---- \n")}

  list(b_hat = sol[1],gamma_hat = sol[2])
}


