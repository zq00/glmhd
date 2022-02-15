#' Setup a system of equations
#'
#' This function sets up a system of four equations for binary regressions.
#' The solution to the system characterizes the asymptotic bias and variance
#' of the M-estimator, as well as the Hessian of the loss function (in case
#' of the MLE, the loss function is the negative log-likelihood).
#'
#' Following is the formula of the four equations:
#' \deqn{
#'  \begin{dcases}
#'  \sigma^2\kappa & =\E{\rho'(S_1)(\lambda\rho'(\mathrm{prox}_{\lambda\rho}(-S_2)))^2+\rho'(-S_1)(\lambda\rho'(\mathrm{prox}_{\lambda\rho}(S_2)))^2}\\
#'  \sigma (1-\kappa) & =\E{\rho'(S_1)Z_2\mathrm{prox}_{\lambda\rho}(\lambda+S_2) + \rho'(-S_1)Z_2\mathrm{prox}_{\lambda\rho}(S_2)}\\
#'  \gamma_0 \alpha & = \E{\rho'(S_1)Z_1\mathrm{prox}_{\lambda\rho}(\lambda+S_2) + \rho'(-S_1)Z_1\mathrm{prox}_{\lambda\rho}(S_2)},
#'  0 & = \E{-\rho'(S_1)\rho'(\mathrm{prox}_{\lambda\rho}(-S_2)) + \rho'(-S_1)\rho'(\mathrm{prox}_{\lambda\rho}(S_2))}.
#'  \end{dcases}
#' }
#' where \eqn{(Z_1, Z_2)\sim\mathcal{N}(0, I_2)} and
#' \deqn{
#'   S_1 = \gamma_0 Z_1 + \beta_0 ,\quad S_2 = \alpha \gamma_0 Z_1 + \sigma Z_2 + b_0,
#' }
#' When the variables does not have an intercept term, then \eqn{b_0 = 0}.
#' If the model does not have an intercept, then \eqn{\beta_0 = 0}.
#'
#' @references
#' \emph{A modern maximum-likelihood theory for high-dimensional logistic regression},
#' Pragya Sur and Emmanuel J. Candes,
#' Proceedings of the National Academy of Sciences Jul 2019, 116 (29) 14516-14525
#' @keywords internal
#' @include prox_op.R integrate2_normal.R
#' @inheritParams find_param
#' @return A function that takes as input the parameters \eqn{(\alpha,\lambda,\sigma,b)}
#'     and returns a vector of length 4, which is the value of the four equations.
#'     When the model contains no intercept term (\code{intercept = FALSE}),
#'     returns a system of three equations. The special case when there is no signal (\code{gamma = 0})
#'     or intercept, returns a system of two equations.
equation_binary <- function(rho_prime, f_prime1, f_prime0, kappa, gamma, beta0, intercept = TRUE){
  function(param, verbose){
    if(gamma == 0) {
      alpha <-  1; lambda <-  param[1]; sigma <-  param[2];
      if(intercept == 0){
        b <-  0
        if(verbose) cat("lambda = ", lambda, ", sigma = ", sigma, "\n")
      }else{
        b <- param[3]
        if(verbose) cat("lambda = ", lambda, ", sigma = ", sigma, ", b = ", b, "\n")
      }
    }else{
      alpha <- param[1]; lambda <-  param[2]; sigma <-  param[3];
      if(intercept == FALSE){
        b <-  0
        if(verbose) cat("alpha = ", alpha, ", lambda = ", lambda, ", sigma = ", sigma, "\n")
      }else{
        b <- param[4]
        if(verbose) cat("alpha = ", alpha, ", lambda = ", lambda, ", sigma = ", sigma, ", b = ", b, "\n")
      }
    }

    s1_fun <- function(x) beta0 + gamma * x[1]
    s2_fun <- function(x) b + alpha * gamma * x[1] + sigma * x[2]

    h1 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      rho_prime(s1) * (s2 - prox_op(f_prime1, lambda, s2))^2 + (1 - rho_prime(s1)) * (s2 - prox_op(f_prime0, lambda, s2))^2
    }

    h2 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      x[2] * prox_op(f_prime1, lambda, s2) * rho_prime(s1) + x[2] * prox_op(f_prime0, lambda, s2) * (1 - rho_prime(s1))
    }

    h3 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      x[1] * prox_op(f_prime1, lambda, s2) * rho_prime(s1) + x[1] * prox_op(f_prime0, lambda, s2) * (1 - rho_prime(s1))
    }

    h4_1 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      f_prime1(prox_op(f_prime1, lambda, s2)) * rho_prime(s1)
    }

    h4_2 <- function(x){
      s1 <- s1_fun(x); s2 <- s2_fun(x)
      f_prime0(prox_op(f_prime0, lambda, s2)) * (1 - rho_prime(s1))
    }

    if(gamma == 0){
      if(intercept == FALSE){
        c(
          sigma^2 * kappa - integrate2_normal(h1),
          sigma* (1 - kappa) - integrate2_normal(h2)
        )
      }else{
        c(
          sigma^2 * kappa - integrate2_normal(h1),
          sigma* (1 - kappa) - integrate2_normal(h2),
          integrate2_normal(h4_1) + integrate2_normal(h4_2)
        )
      }
    }else if(intercept == FALSE){
      c(
        sigma^2 * kappa - integrate2_normal(h1),
        sigma* (1 - kappa) - integrate2_normal(h2),
        gamma * alpha - integrate2_normal(h3)
      )
    }else{
      c(
        sigma^2 * kappa - integrate2_normal(h1),
        sigma* (1 - kappa) - integrate2_normal(h2),
        gamma * alpha - integrate2_normal(h3),
        integrate2_normal(h4_1) + integrate2_normal(h4_2)
      )
    }
  }
}

