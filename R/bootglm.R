#' Compute resized bootstrap MLE for a GLM multiple times
#' 
#' This function implements a parametric bootstrap (for a GLM). Unlike the standard parametric bootstrap,
#' where we use the MLE to generate new responses, here user can input a coefficient vector beta. 
#' The function calls the the function boot_fun, which implements the resized bootstrap once, and returns a vector of MLE.
#' 
#' @param X A covariate matrix of size n * p.
#' @param beta A vector of length n. Coefficients where the parametric bootstrap is applied.
#' @param family An object of class family 
#' @param b_boot An integer of number of bootstrap samples
#' @param verbose Print progress if T
#' 
#' @return  mle_boot a matrix of size p * b_boot of the bootstrap MLE
#' Returns error if the bootstrap MLE does not exist more than 20% of times.
#' @export
bootglm <- function(X, beta, family, b_boot, verbose){
  family$simulate_fun <- glmhd::get_simulate_fun(family) # a function to simulate Y from the linear predictor

  p <- ncol(X)
  b <- 0; nerror <- 0 # b counts bootstrap samples; nerror counts number of times the MLE does not exist
  mle_boot <- matrix(0, p, b_boot) # stores the bootstrap samples
  if(!is.numeric(b_boot)) stop("The number of bootstrap samples must be an integer!")
  if(verbose) cat("\n Total number of bootstrap samples: ", floor(b_boot), "\n")
  while(b < b_boot){
    if(nerror / b_boot > 0.2){error("Bootstrap MLE does not exist for shrinked coefficients!")}
    mle_boot_new <- boot_fun(X, beta, family, sloe = F)
    if(is.numeric(mle_boot_new) && mle_boot_new == -1) {nerror <- nerror + 1; next;}
    b <- b + 1; mle_boot[ ,b] <- mle_boot_new$coef
    
    if(verbose){ if(b %% (b_boot/100) == 0) cat(b, " "); if(b %% (b_boot/10) == 0) cat("\n")}
  }
  
  if(verbose) cat("\n Finished ", b, "resized bootstrap samples. \n")
  
  return(mle_boot)
}
