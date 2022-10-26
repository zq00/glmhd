# bootstrap routine
# generates one parametric bootstrap sample according to the coefficients and glm family given
#     and computes the MLE
# Input: X - covariate matrix 
#     beta - coefficients
#     family - family object for the glm
#     sloe - TRUE if use SLOE to estimate eta
# Output: MLE coefficients
boot_fun <- function(X, beta, family, sloe = F){
  Y <- family$simulate_fun(X%*%beta)
  fit <- tryCatch(error = function(e) {cat("E! "); return(-1)},
                  glm(Y ~ X + 0, family = family))
  if(length(fit) == 1 | !fit$converged){return(-1)}else{
    if(!sloe){
      return(list(coef = fit$coef))
    }else{
      return(list(coef = fit$coef,
                  eta_hat = estimate_eta(X, Y, fit$coef, family))
      )
    }
  }
}
