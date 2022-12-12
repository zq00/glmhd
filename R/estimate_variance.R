#' Estimate std.dev. of the linear predictor evaluated at the MLE
#' 
#' Estimates \eqn{\eta = sd(X_\mathrm{new}^\top \hat{\beta})} where \eqn{X_{\mathrm{new}}} is a new obs. and \eqn{\hat{\beta}} is the MLE
#' when the true coefficient is \eqn{\beta}
#' 
#' Generate \code{b_var} parametric bootstrap samples by sampling new responses \eqn{y} at the
#' observed covariates and when the model coef. is \code{\beta}. 
#' Use SLOE to estimate \eqn{\eta} in each bootstrap sample. 
#' 
#' @inheritParams estimate_eta
#' @param beta A vector of true model coef. 
#' @param b_var Numeric. Number of bootstrap replicates used to estimate \eqn{sd(X^\top \hat{\beta})}.
#' @return A vector of length \code{b_var} of \eqn{\hat{\eta}} in each bootstrap sample. 
#'   if the GLM function reports an error in more than 50% of times, then return -1. 
#'   
#' @export
estimate_variance <- function(x, y, beta, family, b_var){
  nerror <- 0; b <- 0 
  eta_hat <- numeric(b_var) # stores the estimated eta_hat in each bootstrap sample
  while(b < b_var){
    beta_hat_new <- boot_fun(x, beta, family, sloe = F) # Use SLOE to estimate eta
    if(is.numeric(beta_hat_new) && beta_hat_new == -1) {
      nerror <- nerror + 1; # error in fitting a GLM
      if(nerror > 0.5 * b_var) {return(-1)}
      next;
    }
    b <- b + 1; 
    eta_hat[b] <- beta_hat_new$eta_hat # If we want to add other options of how to estimate eta, we can edit here
  }
  
  return(eta_hat)
}
