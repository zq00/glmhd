#' Porportion of linearly separable subsamples
#'
#' This function randomly generate \code{B} subsamples of size \code{nn} and
#' returns the proportion of times that a subsample is linearly separable.
#'
#' @include is_separable.R
#' @family is_separable
#' @param X Covariate matrix. Each row in \code{X} is one observation.
#' @param Y Response vector of \eqn{+1} and \eqn{-1} representing
#'     the two classes. \code{Y} has the same length as the number of rows
#'     in \code{X}.
#' @param B Numeric. Number of subsamples.
#' @param nn Number of observations in each subsample.
#' @return Numeric. The proportion of separable subsamples.
#' @examples
#' n <- 1000; p <- 400
#' X <- matrix(rnorm(n*p, 0, 1), n, p)
#' Y <- 2 * rbinom(n, 1, 0.5) - 1
#' separable_proportion(X, Y, nn = 600, B = 10)
separable_proportion <- function(X, Y, nn, B = 10){
  n <- nrow(X)
  is_sep <- numeric(B)
  for(b in 1:B){
    sind <- sample(1:n, nn, replace = FALSE)
    Xs <- X[sind, ]; Ys = Y[sind]
    is_sep[b] <- is_separable(Xs, Ys)
  }
  mean(is_sep)
}
