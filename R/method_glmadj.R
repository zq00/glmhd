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
#' @export
summary.glmadj <- function(object, ...){
  se <- object$std_adj
  tval <- object$coef_adj / se
  TAB <- cbind(Estimate = object$coef_adj,
               StdErr = se,
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
