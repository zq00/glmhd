#' Estimate eta from MLE coef.
#' 
#' Use an extension of the SLOE estimator to estimate \eqn{\eta = \var(X_{\mathrm{new}}^\top \hat{\beta})^{1/2}},
#' where \eqn{X_{\mathrm{new}}} is a new obs. and \eqn{\hat{\beta}} is the MLE.
#' 
#' Let \eqn{f_y(t)} be the negative log-likelihood when the response is \eqn{y} and linear predictor is \eqn{t}.
#' Define \eqn{w_i = x_i^\top H^{-1}x_i} and \eqn{t_i \ x_i^\top \hat{\beta}}, where \eqn{x_i} is the i-th obs., and
#' \eqn{H} is the Hessian of the negative log-likelihood evaluated at the MLE. Let
#' \deqn{
#' S_i = x_i^\top \hat{\beta} + q_i f'_{y_i}(t_i),
#' }
#' where 
#' \deqn{
#' q_i = \frac{w_i}{1-w_i f''_{y_i}(t_i)}.
#' }
#' Then, the SLOE estimator is defined to be 
#' \deqn{
#' \hat{\eta}^2  = \frac{1}{n}\sum_{i=1}^n S_i^2 - \left(\frac{1}{n}\sum_{i=1}^n S_i\right)^2.
#' }
#'
#'@param X A covariate matrix of size n*p
#'@param y A vector of responses of length n
#'@param beta_hat A vector of MLE (length p)
#'@param family A GLM family, with family and link. See also [compute_deriv()].
#'
#'@return A numeric value of the estimated \eqn{\hat{\eta}}.
#'@export
estimate_eta <- function(X, y, beta_hat, family){
  if(family$family == "binomial") y <- 2 * y - 1
  f <- getg(family)
  g <- f$g
  gprime <- f$gprime
  
  eta_hat <- X %*% beta_hat
  D <- as.vector(gprime(y, eta_hat))
  H <-  t(X) %*% (X * D)
  w <- diag(X %*% (solve(H) %*% t(X)))
  q <- w / (1 - D * w)
  
  eta_tilde <- eta_hat + q * g(y, eta_hat)
  sqrt(mean(eta_tilde^2) - mean(eta_tilde)^2)
}
