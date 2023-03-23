#' Simulate response given linear predictors
#' 
#' @param family A family object, it contains a GLM family \code{family} and an inverse link function \code{linkinv}.
#'   Currently supports binomial and Poisson families. 
#' 
#' @return simulate_fun a function which takes as input a covariate matrix X and a coefficient vector beta.
#' and returns a vector containing simulated responses given X andbeta.
#' 
#' @export
#' @examples
#' family <-  binomial(link = "logit")
#' fun <- get_simulate_fun(family)
#' y <- fun(X = matrix(rnorm(100, 0, 1), 10, 10), beta = rnorm(10))
get_simulate_fun <- function(family){
  if(family$family == "binomial"){
    simulate_fun <- function(X, beta){
      t <- X %*% beta
      stats::rbinom(length(t), 1, family$linkinv(t))
    }
  }
  if(family$family == "poisson"){
    simulate_fun <- function(X, beta){
      t <- X %*% beta
      stats::rpois(length(t), family$linkinv(t))
    }
  }
  return(simulate_fun)
}