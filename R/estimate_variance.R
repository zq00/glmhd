#' Estimate std.dev. of the linear predictor evaluated at the MLE
#' 
#' Estimates \eqn{\eta = sd(X_\mathrm{new}^\top \hat{\beta})} where \eqn{X_{\mathrm{new}}} is a new obs. and \eqn{\hat{\beta}} is the MLE
#' when the true coefficient is \eqn{\beta}
#' 
#' Generate \code{b_var} parametric bootstrap samples by sampling new responses \eqn{y} at the
#' observed covariates and when the model coef. is \code{\beta}. 
#' Use SLOE to estimate \eqn{\eta} in each bootstrap sample. 
#' 
#' @param x A covariate matrix of size n*p.
#' @param beta A vector of true model coef. (used to sample \code{y})
#' @param family A GLM family, with family and link. See also [compute_deriv()].
#' @param b_var Numeric. Number of bootstrap samples
#' @return A vector of length \code{b_var} of \eqn{\hat{\eta}} in each bootstrap sample. 
#'   if the GLM function reports an error in more than 50% of times, then return -1 (which means the MLE may not exist at \eqn{\beta}). 
#'   
#' @export
estimate_variance <- function(x, beta, family, b_var){
  nerror <- 0; b <- 0 
  eta_hat <- numeric(b_var) # stores the estimated eta_hat in each bootstrap sample
  while(b < b_var){
    beta_hat_new <- boot_fun(x, beta, family,sloe = T) # Use SLOE to estimate eta
    if(is.numeric(beta_hat_new) && beta_hat_new == -1) {
      nerror <- nerror + 1; # error in fitting a GLM
      if(nerror > 0.5 * b_var) {return(-1)}
      next;
    }
    b <- b + 1; 
    eta_hat[b] <- beta_hat_new$eta_hat 
  }
  
  return(eta_hat)
}
