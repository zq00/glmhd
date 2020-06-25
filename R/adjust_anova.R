#' LRT for high-dimensional glm
#'
#' This is the generic function to compute the likelihood ratio statistics
#' and p-valus for binary regression models.
#'
#' Under the null hypothesis:
#' \deqn{
#' H_0:\quad \beta_1 = \ldots = \beta_k = 0,
#' }
#' two times the likelihood ratio statistics \eqn{\Lambda} is asymptotically from a
#' re-scaled chi-squared distribution with \eqn{k} degrees of freedom
#' \deqn{
#' 2 \Lambda \stackrel{d}{\longrightarrow} \frac{\sigma_\star^2}{\lambda_\star}\chi^2_k
#' }
#' where \eqn{(\sigma_\star, \lambda_\star)} can be computed using \code{find_param} function.
#'
#' @references
#' \emph{The likelihood ratio test in high-dimensional logistic regression is asymptotically a rescaled Chi-square},
#' Pragya Sur, Yuxin Chen and Emmanuel Candes, Probab. Theory Relat. Fields 175, 487–558 (2019).
#'
#' @param object A list with at least two elements from the class \link[stats]{glm}.
#'     The function does not check whether the models are nested.
#' @param param A named vector with elements \code{sigma_s} and \code{lambda_s} which
#'     can be outputs from the \code{adjust_glm} function.
#' @return A list with two elements
#' \describe{
#' \item{models}{Model formula}
#' \item{anova.tab}{Significance tables for the LRT}
#' }
#' @importFrom stats pchisq
#' @examples
#' n <- 1000L
#' p <- 300L
#' X <- matrix(rnorm(n*p, 0, 1), n, p) / sqrt(p)
#' beta <- rep(c(0, 1), each = p / 2)
#' Y <- rbinom(n, 1, 1 / (1 + exp(- X %*% beta)))
#' f1 <- glm(Y ~ X + 0, family = binomial, x = TRUE, y = TRUE)
#' f2 <- glm(Y ~ X[ , -1] + 0, family = binomial, x = TRUE, y = TRUE)
#' adjusted_fit <- adjust_glm(f1)
#' lrt_glm(list(f1, f2), adjusted_fit$param)
#' @export
lrt_glm <- function(object, param){
  if(length(object) < 2){
    stop("Please input at least two glm models!")
  }else{
    # Extract the log-likelihood and degree of freedom
    df <- sapply(object, function(x) x$df.residual)
    dev <- sapply(object, function(x) x$deviance)
    # Order models from small to large
    ord <- sort(df, decreasing = TRUE, index.return = TRUE)$ix
    # Compute the deviance
    diff_df <- - diff(df[ord])
    diff_dev <- - diff(dev[ord])
    # Compute the p-value
    p_val <- pchisq(diff_dev * param["lambda_s"] / param["sigma_s"]^2,
                    df = diff_df, lower.tail = FALSE)
    # Significance table
    tab <- data.frame("Resid Df" = df[ord],
                 "Resid Dev" = dev[ord],
                 "Df" = c(NA, diff_df),
                 "Deviance" = c(NA, diff_dev),
                 "p.value" = c(NA, p_val))
    rownames(tab) <- 1:length(object)
    # Models
    mod <- list()
    for(i in 1:length(object)){
      mod[[i]] <- object[[i]]$formula
    }
  }
  result <- list(
    models = mod,
    anova.tab = tab
  )
  class(result) <- "lrt_adj"
  result
}

print.lrt_adj <- function(object){
  cat("Analysis of Deviance Table \n")
  for(i in 1:length(object$models)){
    cat(paste0("Model ", i, ": "))
    print(object$models[[i]])
  }
  printCoefmat(object$anova.tab, P.values = TRUE, has.Pvalue = TRUE, na.print = " ")
}
