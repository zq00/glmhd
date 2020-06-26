#' Methods for adjusted glm objects
#'
#' Availabel methods for \code{glmadj} objects: print coefficients, summary,
#' prediction
#' @name glmadj
#' @param x A \code{glmadj} object, created from \code{\link{adjust_glm}}.
#' @param ... Additional parameters
#' @importFrom stats model.matrix printCoefmat
#' @seealso \code{\link{adjust_glm}}
#' @export

print.glmadj <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coef_adj)
}

#' @rdname glmadj
#' @param object A \code{glmadj} object, created from \code{\link{adjust_glm}}.
#' @return A summary table with four columns
#' \itemize{
#' \item ajusted_mle \eqn{\hat{\beta}_j^\mathrm{Adj} = \hat{\beta}_j^\mathrm{MLE} / \alpha_\star}.
#' \item adjusted_std. Standard error of the adjusted MLE \eqn{\hat{\sigma}_j / \alpha_\star = \sigma_\star / \alpha_\star \tau_j}.
#' \item t.value \eqn{t_j = \hat{\beta}_j^\mathrm{MLE} / \hat{\sigma}_j}. When \eqn{\beta_j = 0}, \eqn{t_j} is approximately
#'     a standard Gaussian as \eqn{n,p \to\infty}.
#' \item p.value 2-sided p-value using t.value to test whether \eqn{\beta_j = 0}.
#' }
#' @export
summary.glmadj <- function(object, ...){
  se <- object$std_adj / object$param["alpha_s"])
  tval <- object$coef_adj / se
  TAB <- cbind(adjusted_mle = object$coef_adj,
               adjusted_std = se,
               t.value = tval,
               p.value = 2 * pnorm(-abs(tval)))
  result <- list(call=object$call,
                 coefficients=TAB)
  class(result) <- "summary.glmadj"
  result
}

#' @rdname glmadj
#' @export
print.summary.glmadj <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE, ...)
}

#' @param newdata A data frame created with \link[stats]{model.frame}.
#' @rdname glmadj
#' @export
predict.glmadj <- function(object, newdata = NULL, ...){
  if(is.null(newdata))
    y <- object$fitted.values
  else{
    if(!is.null(object$glm_output$formula)){
      # model has been fitted using formula interface
      x <- model.matrix(object$glm_output$formula, newdata)
    }
    else{
      warning("Missing model formula!")
      x <- newdata
    }
    if(!is.null(object$intercept)){
      y <- as.vector(object$intercept + x[ ,-1] %*% object$coef_adj)
    }else{
      y <- as.vector(x %*% object$coef_adj)
    }
  }
  y
}
