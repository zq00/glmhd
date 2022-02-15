#' Compute the phase transition curve
#'
#' \code{solve_kappa} computes the problem dimension \eqn{\kappa} where
#' the phase transition occurs in binary regression, given \eqn{\beta} and \eqn{\gamma_0}. \code{solve_beta} and \code{solve_gamma}
#' computes \eqn{\beta_0} and \eqn{\gamma_0} on the phase transition curve
#' given the other one and \eqn{\kappa}.
#'
#' When covariates are multivariate Gaussian, the phase transition dimension
#' can be characterized as following.
#' \deqn{
#' \kappa > h_{\mathrm{MLE}}(\beta_0, \gamma_0) \implies \lim_{n,p\to\infty} \mathrm{P}(\text{MLE exists}) = 0
#' }
#' \deqn{
#' \kappa < h_{\mathrm{MLE}}(\beta_0, \gamma_0) \implies \lim_{n,p\to\infty} \mathrm{P}(\text{MLE exists}) = 1.
#' }
#' The function \eqn{h} is defined to be
#' \deqn{
#' h_{\mathrm{MLE}}(\beta_0, \gamma_0) = \min_{t_0, t_1 \in \mathbb{R}} \mathbb{E}\left[(t_0 Y + t_1 V - Z)_+^2 \right],
#' }
#' where \eqn{X\sim\mathcal{N}(0,1)} and \eqn{\mathrm{P}(Y=1|X) = 1- \mathrm{P}(Y=-1|X) = \rho'(\beta_0 + \gamma_0 X) }.
#' \eqn{Z\sim\mathcal{N}(0,1)} and is independent of \eqn{X,Y}. The phase transition
#' curve is thus \eqn{\kappa(\beta_0, \gamma_0)}. It also depends on the
#' success probability \eqn{\rho'}.
#'
#' @include hinge.R integrate2_normal.R
#' @importFrom stats uniroot optim
#' @param rho_prime Function. Success probability \eqn{\rho'(t) = \mathrm{P}(Y=1\,|\, X^\top \beta = t)}
#' @param beta0 Numeric. Intercept value.
#' @param gamma0 Numeric. Signal strength.
#' @param kappa Numeric. Problem dimension on the phase transition curve.
#' @param verbose Print progress if \code{TRUE}.
#' @return Numeric. Problem dimension \eqn{\kappa} (\eqn{\beta} or \eqn{\gamma}) on the phase transition curve.
#' @references
#' \emph{The phase transition for the existence of the maximum likelihood estimate in high-dimensional logistic regression}
#' Emmanuel J. Candes and Pragya Sur, Ann. Statist., Volume 48, Number 1 (2020), 27-42.
#'
#' @examples
#' \dontrun{
#' # when Y is independent of X, should return 0.5 for logistic model
#' # should return 0.5
#' rho_prime_logistic <- function(t) 1 / (1 + exp(-t))
#' solve_kappa(rho_prime_logistic, 0, 0)
#' }

solve_kappa <- function(rho_prime, beta0, gamma0){
  h <- function(t){
    f1 <- function(x) (hinge(t[1] + t[2] * x[1] - x[2]))^2 * rho_prime(beta0 + gamma0 * x[1])
    f2 <- function(x) (hinge(-t[1] - t[2] * x[1] - x[2]))^2 * (1 - rho_prime(beta0 + gamma0 * x[1]))
    integrate2_normal(f1) + integrate2_normal(f2)
  }

  opt <- optim(par=c(0,0), h, method = "BFGS", control = list(abstol = 1e-5, maxit = 200))
  return(list(conv = opt$convergence, val = opt$val))
}

#' @rdname solve_kappa
solve_beta <- function(rho_prime, kappa, gamma0, verbose = FALSE){
  if(solve_kappa(rho_prime, 0, gamma0) < kappa){
    return(-1)
  }else{
    f <- function(beta){
      val <- solve_kappa(rho_prime, beta, gamma0) - kappa
      if(verbose) cat("beta = ", beta, "; diff = ", val, "\n")
      val
    }
    beta_hat <- uniroot(f, interval = c(0, 10), extendInt = "yes", tol = 0.001)$root
    if(verbose) cat("beta_hat = ", beta_hat, "\n")
    return(beta_hat)
  }
}

#' @rdname solve_kappa
solve_gamma <- function(rho_prime, kappa, beta0, verbose = FALSE){
  if(solve_kappa(rho_prime, beta0, 0) < kappa){
    return(0)
  }else{
    f <- function(gamma) {
      val <- solve_kappa(rho_prime, beta0, gamma) - kappa
      if(verbose) cat("gamma = ", gamma, "; diff = ", val, "\n")
      val
    }
    gamma_hat <- uniroot(f, interval = c(0, 10), extendInt = "yes", tol = 0.001)$root
    if(verbose) cat("gamma_hat = ", gamma_hat, "\n")
    return(gamma_hat)
  }
}

