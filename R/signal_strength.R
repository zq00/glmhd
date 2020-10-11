#' Estimate intercept and signal strength
#'
#' \code{signal_strength} estimates \eqn{(\beta_0, \gamma)} using the estimated \eqn{\kappa_s}
#' and observed proportion of successes.
#'
#' Assume that \eqn{Y} depends on \eqn{X} as
#' \deqn{
#' \mathrm{P}(Y=1\,|\,X) = \rho'(X^\top \beta + \beta_0),
#' }
#' and let the signal strength be \eqn{\gamma = \mathrm{Var}(X^\top \beta)^{1/2}}.
#' The pair \eqn{(\beta_0, \gamma)} satisfies
#' \itemize{
#' \item They are on the phase transition curve \eqn{\kappa(\beta_0, \gamma)}.
#'     \deqn{
#'     \hat{\kappa} \approx \kappa(\beta_0, \gamma)
#'      }
#' \item The observed proportion of \eqn{Y=1} should be close to the expected proportion.
#'     \deqn{
#'     p_0 \approx \prob(Y = 1\,|\, \beta_0, \gamma) = \mathrm{E}{\mathrm{Ber}(\rho'(\beta_0 + \gamma Z)}
#'     }
#'    where \eqn{Z} is a standard normal variable.
#' }
#' We solve the above system of two equations to obtain an estimate of \eqn{(\beta_0, \gamma)}
#'
#' @importFrom pracma fsolve
#' @importFrom stats integrate
#' @include prox_op.R h_eq.R find_param.R
#' @param rho_prime A function that computes the success probability \eqn{\rho'(t) = \mathrm{P}(Y=1 | X^\top \beta = t)},
#'     here \eqn{\beta} is the coefficient. The default is logistic model.
#' @param kappa_hat Numeric. Estimated dimension where the data becomes linearly separable
#' @param intercept Logical \code{TRUE} if the model contains an intercept.
#' @param p0 Numeric. Proportion of outcomes \eqn{Y=1}.
#' @param verbose Should progress be printed? If \code{TRUE}, prints progress at each step.
#' @param tol Numeric. Tolerance to be used in \code{fsolve} function.
#'
#' @return If the model does not contain an intercept, returns estimated \code{gamma_hat}.
#'     Otherwise, returns a list with two components
#' \describe{
#'   \item{gamma_hat}{Estimated signal strength.}
#'   \item{b_hat}{Estimated intercept.}
#' }
#' @examples
#' \dontrun{
#' # no signal case
#' # should return 0, returns 0.0127
#' signal <- signal_strength(kappa_hat = 0.5, intercept = FALSE)
#' signal$gamma_hat
#' }
signal_strength <- function(rho_prime = rho_prime_logistic,
                            kappa_hat,
                            intercept = FALSE,
                            p0 = NA,
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


