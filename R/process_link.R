#' Extract likelihood from link function
#'
#' This function returns success probability and derivatives of the log-likelihood
#' for commonly used link function in the GLM binary family.
#'
#' @keywords internal
#' @param glm_family A glm \code{familly} object for \code{binomial} family.
#'     It should have an element called \code{link}.
#' @return A list with three elements
#' \describe{
#' \item{rho_prime}{A function that computes the success probabily \eqn{\rho'(t) = \mathrm{P}(Y=1|X^\top \beta = t)}.}
#' \item{f_prime1}{Derivative of the negative log-likelihood when \eqn{Y=1}}
#' \item{f_prime0}{Derivative of the negative log-likelihood when \eqn{Y=0}}
#' }
#' @importFrom stats pnorm pcauchy dnorm dcauchy
#' @examples
#' \dontrun{
#' process_link(binomial(link = "cloglog"))
#' }
process_link <- function(glm_family){
  rho_prime <-  switch(glm_family$link,
                       "logit" = rho_prime_logistic,
                       "probit" = function(t) pnorm(t),
                       "cauchit" = function(t) pcauchy(t),
                       "log" = function(t) exp(t),
                       "cloglog" = function(t) 1 - exp(-exp(t))
  )
  f_prime1 <- switch(glm_family$link,
                     "logit" = f_prime1_logistic,
                     "probit" = function(t) - dnorm(t) / pnorm(t),
                     "cauchit" = function(t) - dcauchy(t) / pcauchy(t),
                     "log" = function(t) -t,
                     "cloglog" = function(t) - exp(t) / (exp(exp(t)) - 1)
  )
  f_prime0 <- switch(glm_family$link,
                     "logit" = f_prime0_logistic,
                     "probit" = function(t) dnorm(t) / pnorm(-t),
                     "cauchit" = function(t) dcauchy(t) / pcauchy(-t),
                     "log" = function(t) 1 / (exp(-t) - 1),
                     "cloglog" = function(t) exp(t)
  )
  list(rho_prime = rho_prime,
       f_prime1 = f_prime1,
       f_prime0 = f_prime0)
}
