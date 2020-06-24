#' Estimate coefficients of a high-dimensional GLM
#'
#' This is the generic function and methods associated with the adjusted MLE for
#' high dimensional generalized linear models (glm).
#'
#' The function \code{adjust_glm} is built on the \code{\link[stats]{glm}} function: the input should be from the
#' class \link[stats]{glm}, such as the output from \code{glm}.
#' Currently, it works with binary regressions only, for which it calls the function \code{adjust_binary}
#' to estimate coefficients and their standard deviations. In high dimensions,
#' the magnitude of MLE is biased upward, and the variance estimates from classical theory based on
#' inverse Fisher information underestimates the true variability. \code{adjust_binary} adjusts the MLE and variance estimates
#' applying the asymptotic theory assuming \eqn{p, n\to\infty} (\eqn{p<n}), and is exact when
#' covariates \eqn{X} is multivariate Gaussian.
#'
#' \code{adjust_glm} creates an object from class \code{glmadj}. \code{summary.glmadj} creates a data table of adjusted coefficients,
#' estiamtes of the standard error and p-values using a two-sided t-test. If
#' the covariates are multivariate Gaussian, the p-values are asymptotically
#' from a uniform distribution.
#'
#' \code{predict.glmadj} computes the prediction \eqn{X^\top \hat{\beta}^{\mathrm{Adj}}} for
#' new input data. The new data should be of the same format which you
#' used to fit the \code{glm}. If there's no new data, it outputs the estimated
#' \eqn{\eqn{X^\top \hat{\beta}_{\mathrm{Adj}}}} for the original dataset.
#'
#' @param glm_output An object from the class \link[stats]{glm},
#'     for example the output from the \code{\link[stats]{glm}} function
#' @param verbose If \code{TRUE}, print progress at every step.
#' @param echo If \code{TRUE}, returns the original \code{glm_output}.
#' @seealso \code{\link{adjust_binary}}
#' @references [TODO: add]
#' @examples
#' # Problem size
#' n <- 1000L
#' p <- 300L
#' # Generate data matrix
#' X <- matrix(rnorm(n*p, 0, 1), n, p)
#' # Generate parameters
#' beta <- rep(c(0, 1), p / 2) / sqrt(p) * 2
#' #' ## No intercept model
#' # Sample response vector
#' Y <- rbinom(n, 1, 1 / (1 + exp(- X %*% beta)))
#' # Fit a glm, by default computes a logistic regression
#' fit <- glm(Y ~ X + 0, family = binomial, x = TRUE, y = TRUE)
#' adjusted_fit <- adjust_glm(fit)
#' # Print summary
#' summary(adjusted_fit)
#' # Predict on new data
#' predict(adjusted_fit)
#' # Extract adjusted coefficients
#' adjusted_fit$coef_adj
#' # Extract adjusted standard error
#' adjusted_fit$std_adj
#' @export
adjust_glm <- function(glm_output, verbose, echo, ...) UseMethod("glmadj")

#' @rdname adjust_glm
glmadj.default <- function(glm_output, verbose = FALSE, echo = TRUE, ...){
  if(glm_output$family$family == "binomial"){
    mle_adj <- adjust_binary(glm_output, verbose, echo)
  }else{
    stop("Not binary regression!")
  }
  # Fitted value is the linear part, i.e. eta
  if(!is.null(mle_adj$intercept)){
    mle_adj$fitted.values <- as.vector(mle_adj$intercept +
                                         glm_output$x[ ,-1] %*% mle_adj$coef_adj)
  }else{
    mle_adj$fitted.values <- as.vector(glm_output$x %*% mle_adj$coef_adj)
  }
  mle_adj$call <- match.call()
  class(mle_adj) <- "glmadj"
  mle_adj
}

#' @rdname adjust_glm
print.glmadj <- function(object, ...){
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(object$coef_adj)
}

#' @rdname adjust_glm
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

#' @rdname adjust_glm
print.summary.glmadj <- function(object, ...){
  cat("Call:\n")
  print(object$call)
  cat("\n")
  printCoefmat(object$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

#' @rdname adjust_glm
predict.glmadj <- function(object, newdata = NULL, ...){
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$glm_output$formula)){
      # model has been fitted using formula interface
      x <- model.matrix(object$glm_output$formula, newdata)
    }
    else{
      warning("Missing model formula!")
      x <- newdata
    }
    if(!is.null(mle_adj$intercept)){
      y <- as.vector(object$intercept + x[ ,-1] %*% object$coef_adj)
    }else{
      y <- as.vector(x %*% object$coef_adj)
    }
  }
  y
}


















