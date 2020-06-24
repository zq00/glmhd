#' Estimate the problem dimension where two classes become linearly separable
#'
#' This function estimates the sample size \eqn{n_s}, or equivalently problem dimension
#' \eqn{\kappa_s = p/n_s}, that two classes from the data becomes separable. To locate \eqn{\kappa_s},
#' we bisect the interval \eqn{[p/n, 0.5]}, until the window
#' size is smaller than \code{eps}. For each sample size \code{nn}, it generates
#' \code{B} subsamples of size \code{nn}, and estimate the separable probability
#' \eqn{\hat{\pi}} with the proportion of separable subsamples.
#' Finally we fit a logistic regression using \eqn{\hat{\pi}} as response
#' and \eqn{\kappa = p/nn} as covariate to determine the \eqn{\hat{\kappa}}
#' where separable probability is 0.5.
#'
#' @include is_separable.R separable_proportion.R
#' @param X Covariate matrix. Each row in \code{X} is one observation.
#' @param Y Response vector of \eqn{+1} and \eqn{-1} representing
#'     the two classes. \code{Y} has the same length as the number of rows
#'     in \code{X}.
#' @param B Numeric. How many subsamples should I generate
#'     for each sample size?
#' @param eps Numeric. Minimum window size.
#'     Terminate when the search interval is smaller than \code{eps}
#' @param verbose Print prgress if \code{TRUE}.
#' @return Numeric. Estimated kappa_hat.
#' @examples
#' # Y is independent of X, kappa_s is approximately 0.5
#' n <- 1000; p <- 200
#' X <- matrix(rnorm(n*p, 0, 1), n, p)
#' Y <- 2 * rbinom(n, 1, 0.5) - 1
#' probe_frontier(X, Y, verbose = TRUE)
probe_frontier <- function(X, Y, B = 10, eps = 0.001, verbose = FALSE){
  # Problem dimension
  n <- nrow(X); p <- ncol(X)
  kappa <- p / n
  # Binary search for kappa_s
  l <- kappa
  u <- 0.51
  k_new <- (l + u) / 2
  pi_hats <- NULL # Proportion of separable samples
  if(verbose) cat("------ Begin Probe Frontior ------\n")
  while(u - l > eps){
    nn = floor(p / k_new) # Sample size in each subsample
    # Estimate separable probability
    pi_new <- separable_proportion(X, Y, nn, B)
    if(verbose) cat("kappa = ", k_new, "; pi_hat = ", pi_new, "\n")
    pi_hats <- rbind(pi_hats, c(k_new, pi_new))
    # Update kappa_hat
    if(pi_new>0.5){
      u <- k_new
    }else{
      l <- k_new
    }
    k_new <- (l+u)/2
  }
  # Estimate kappa_s
  fit <- glm(pi_hats[ ,2] ~ pi_hats[ ,1], family = binomial, weights = rep(B,nrow(pi_hats)))
  kappa_s <-  - fit$coef[1] / fit$coef[2]
  if(verbose) cat("Found kappa_s = ", ifelse(kappa_s > 0.5, 0.5, kappa_s), "\n")
  if(verbose) cat("------- \n")
  if(kappa_s > 0.5) {
    return(0.5)
  }else{
    return(kappa_s)
  }
}
