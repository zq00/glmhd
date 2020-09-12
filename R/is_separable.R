#' Are the two classes linearly separable?
#'
#' The function \code{is_separable} determines if the two classes \eqn{Y=1}
#' and \eqn{Y=-1} are linearly separable using columns of \eqn{X}.
#'
#' The two classes are separable if there exists a non-zero vector \eqn{b}
#' that satisfies \deqn{y_i x_i^\top b \leq 0} for every observation \eqn{i}.
#'
#' @param X Covariate matrix. Each row of \code{X} corresponds to one observation.
#' @param Y A vector of \eqn{+1} and \eqn{-1} representing the two classes.
#'     \code{Y} has the same length as the number of rows in \code{X}.
#' @param add_intercept If \code{TRUE}, add an intercept to the matrix \code{X}.
#'     Set to \code{FALSE} by default.
#' @return Returns 1 if the data is separable and 0 otherwise.
#'
#' @keywords internal
#' @importFrom ECOSolveR ECOS_csolve
#' @importFrom Matrix rankMarix
#' @export
#' @seealso \link[safeBinaryRegression]{glm}
#' @references
#' \emph{The phase transition for the existence of the maximum likelihood estimate in high-dimensional logistic regression}
#' Emmanuel J. Candes and Pragya Sur, Ann. Statist., Volume 48, Number 1 (2020), 27-42.
#'
#' \emph{Linear programming algorithms for detecting separated data in binary logistic regression models}, Kjell Konis, Ph.D. thesis, Univ. Oxford.
#' @examples
#' \dontrun{
#' n <- 1000; p <- 400
#' X <- matrix(rnorm(n*p, 0, 1), n, p)
#' Y <- 2 * rbinom(n, 1, 0.5) - 1
#' is_separable(X, Y, add_intercept = TRUE)
#' }
is_separable <- function(X, Y, add_intercept = FALSE){
  # Problem dimension
  n <- nrow(X); p <- ncol(X)

  if(add_intercept){
    Z <- - cbind(1, X) * Y
  }else{
    Z <- - X * Y
  }

  if(rankMatrix(Z, method = "qrLINPACK")[1] < p) return(1)

  zbar <- colMeans(Z)
  # Is 0 solution to the LP?
  c <-rep(0, n)
  G <- diag(n)
  h <- rep(0, n)
  A <- t(Z)
  b <- zbar
  dims <- list(l = n, q = NULL, e = 0L) # l is dimension of positive orthant cones
  # solves minimize c^t x, such that Ax = b, h - Gx in K
  sol <- ECOS_csolve(c = c, G = G, h = h, dims = dims, A = A, b = b)

  # If infeasible, then the MLE does not exist
  ifelse(sol$summary["pinf"] == 1, return(1), return(0))
}
