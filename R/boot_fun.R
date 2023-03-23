#' Compute one bootstrap for a GLM
#' 
#' Generates one parametric bootstrap sample given GLM family and model coefficients.
#' 
#' This model does not fit an intercept. As a result, if you would like to fit an intercept,
#' please include a column of 1 in the matrix X. 
#' 
#' @param X A n*p covariate matrix (each row represents one obs.).
#' @param beta Model coefficient.
#' @param family GLM family. Besides the usual definition of glm family (see [stats::family]), 
#'   \code{family} should include an additional element \code{simulate_fun} (see [get_simulate_fun]).
#' @param sloe If True, use SLOE to estimate \eqn{\eta = \var(X_{\mathrm{new}^\top \hat{\beta}})^{1/2}} 
#'   (see the function \code{estimate_eta}).
#' 
#' @return 
#' \describe{
#' \item{coef}{A length p vector of MLE coefficients.}
#' \item{eta_hat}{If \code{sloe = T}, returns the estimated \eqn{\hat{\eta}}.}
#' }
#' @export

boot_fun <- function(X, beta, family, sloe = F){
  
  Y <- family$simulate_fun(X, beta)
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
