#' Simulate response given linear predictors
#' 
#' @param family A family object, it contains a GLM family \code{family} and an inverse link function \code{linkinv}.
#'   Currently supports binomial and Poisson families. 
#' 
#' @return simulate_fun a function which takes a vector of linear predictors as inputs
#' and returns a vector of simulated responses given the linear predictors.
#' 
#' @export
#' @examples
#' family <-  binomial(link = "logit")
#' fun <- get_simulate_fun(family)
#' y <- fun(rnorm(10))
get_simulate_fun <- function(family){
  if(family$family == "binomial"){
    simulate_fun <- function(t) stats::rbinom(length(t), 1, family$linkinv(t))
  }
  if(family$family == "poisson"){
    simulate_fun <- function(t) stats::rpois(length(t), family$linkinv(t))
  }
  return(simulate_fun)
}