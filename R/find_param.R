#' Solves a system of nonlinear equations
#'
#' This function solves a system of equations, whose solution characterizes
#' the asymptotic bias and variance of the M-estimator (in case of the MLE, it is the negative log-likelihood).
#'
#' @param rho_prime A function that computes the success probability \eqn{\rho'(t) = \mathrm{P}(Y=1 | X^\top \beta = t)},
#'     here \eqn{\beta} is the coefficient. The default is logistic model.
#' @param f_prime1 A function. Derivative of the loss function when \eqn{Y = 1}.
#'     The default is the derivative of the negative log-likelihood of logistic
#'     regression when \eqn{Y = 1}.
#' @param f_prime0 A function. Derivative of the loss function when \eqn{Y = -1}.
#'     The default is the derivative of the negative log-likelihood of logistic
#'     regression when \eqn{Y = -1}.
#' @param kappa Numeric. The problem dimension \eqn{\kappa = p/n}.
#' @param gamma Numeric. Signal strength \eqn{\gamma = \sqrt{\mathrm{Var}(X^\top \beta)}}.
#' @param beta0 Numeric. Intercept.
#' @param intercept If \code{TRUE}, the glm contains an intercept term.
#'     \code{intercept = TRUE} by default.
#' @param verbose If \code{TRUE}, print progress at each step.
#' @return A vector solution to the system. When \code{gamma != 0} and \code{b !=0},
#'     returns \eqn{(\alpha_\star, \lambda_\star, \sqrt{\kappa}\sigma_\star, b_\star)}.
#'     When signal strength is zero (\code{gamma = 0}), returns the solution to the system with
#'     three equations \eqn{(\lambda_\star, \sqrt{\kappa}\sigma_\star, b_\star)}. When \code{gamma = 0} and
#'     \code{b = 0}, returns \eqn{(\lambda_\star, \sqrt{\kappa}\sigma_\star)}.
#' @importFrom pracma fsolve
#' @include equation_binary.R prox_op.R integrate2_normal.R logistic_model.R
#' @references
#' \emph{The Impact of Regularization on High-dimensional Logistic Regression},
#' Fariborz Salehi, Ehsan Abbasi and Babak Hassibi, Proceedings of NeurIPS 2019.
#' @export
#' @examples
#' # Compute parameters for a logistic model
#' param <- find_param(kappa = 0.1, gamma = sqrt(5))
#' # Asymptotic bias
#' param[1]
#' # Standard deviation
#' param[3] / sqrt(0.1)
#' # Another example
#' param <- find_param(kappa = 0.1, gamma = 0, intercept = FALSE)
#' # Asymptotic standard deviation
#' param[2] / sqrt(0.1)

find_param <- function(rho_prime = rho_prime_logistic,
                       f_prime1 = f_prime1_logistic,
                       f_prime0 = f_prime0_logistic,
                       kappa,
                       gamma,
                       beta0 = 0,
                       intercept = TRUE,
                       verbose = FALSE){
  if(gamma == 0){ # case of no signal
    if(intercept == FALSE){
      x_init <- c(2, 2)
    }else{
      x_init <- c(2, 2, beta0)
    }
  }else if(intercept == FALSE){
    x_init <- c(2, 2, 1 + gamma * 2)
  }else{
    x_init <- c(2, 2, 1 + gamma * 2, beta0)
  }
  # Setup system of equations
  f_eq <- equation_binary(rho_prime, f_prime1, f_prime0,
                          kappa, gamma, beta0, intercept)
  if(verbose) {
    cat("Solve parameters for: kappa = ", kappa, ", gamma = ", gamma, ", beta0 = ", beta0, "\n")
    if(intercept == TRUE) cat(", with intercept. \n")
  }
  sol <- fsolve(f_eq, x_init, J = NULL, maxiter = 100, tol = 1e-4, verbose)
  if(verbose) cat("Solution is", sol$x)

  if(gamma == 0){c(1,sol$x)}else{sol$x}
}

