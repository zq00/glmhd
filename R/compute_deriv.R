#' Derivative of the negative log-likelihood 
#' 
#' @description 
#' Compute the derivative of the negative log-likelihood of a GLM w.r.t. the linear predictor
#' 
#' Let \eqn{f_y(t)} be the negative log-likelihood function when the linear predictor is \eqn{t}
#' and the response is \eqn{y}. This function computes 
#' \deqn{
#' g(t) = f'_y(t)
#' }
#' and 
#' \deqn{
#' g'(t)=f''_y(t).
#' }
#' This function uses *formula* of GLM likelihoods. Currently, it supports Poisson regression (log link) and binary regression 
#' (with logit/probit link).
#' 
#' For a logistic regression, \eqn{f_y(t) = \log(1+e^{-yt})}, where \eqn{y \in \pm 1}. Then
#' \deqn{
#' g(t) = \frac{-y}{1+e^{yt}},
#' }
#' and 
#' \deqn{
#' g'(t) = \frac{1}{(1+e^{yt})(1+e^{-yt})}. 
#' }
#' 
#' For a probit regression, \eqn{f_y(t) = -\log(\Phi(yt))} where \eqn{\Phi(\cdot)} is the normal cdf.
#' Then,
#' \deqn{
#' g(t) = -\frac{y\phi(yt)}{\Phi(yt)},
#' }
#' where \eqn{\phi} is the normal pdf, and 
#' \deqn{
#' g'(t) = \frac{\phi(yt)^2}{\Phi(yt)^2} - \frac{\phi'(yt)}{\Phi(yt)}.
#' }
#' For a Poisson regression, \eqn{f_y(t) = e^t - yt + \log(y!)}, so 
#' \deqn{
#' g(t) = e^t - y,
#' }
#' and 
#' \deqn{
#' g't(t) = e^t. 
#' }
#' 
#' @param family A GLM family with family and link 
#' @return A list of two functions
#' \describe{
#' \item{g}{A function which takes two inputs -- a vector of response y and a vector of linear predictors t -- and returns the derivative of the negative log-likelihood  w.r.t. the linear predictor at y}
#' \item{gprime}{A function which takes two inputs -- a vector of response y and a vector of linear predictors t -- and returns the *second* derivative of the negative log-likelihood  w.r.t. the linear predictor at y}
#' }
#' @export
compute_deriv <- function(family){
  if(family$family == "poisson"){
    if(family$link == "log"){
      g <- function(y, t) -y + exp(t)
      gprime <- function(y, t) exp(t)
    }
  }
  
  if(family$family == "binomial"){
    if(family$link == "logit"){
      g <- function(y, t) -y/(1 + exp(y * t))
      gprime <- function(y, t) 1/(1+exp(t))/(1+exp(-t))
    }else if(family$link == "probit"){
      g <- function(y, t) -y * dnorm(y*t)/pnorm(y * t)
      gprime <- function(y, t) {
        phiprime <- -(y * t)*exp(-t^2 / 2)/sqrt(2*pi)
        dnorm(y*t)^2 / pnorm(y * t)^2 - phiprime / pnorm(y*t)
      }
    }else if(family$link == "cloglog"){ # not verified
      g <- function(y, t){
        n <- length(y)
        val <- numeric(n)
        for(i in 1:n){
          val[i] <-ifelse(y[i] > 0, - exp(t[i]) / (exp(exp(t[i])) - 1), exp(t[i]))
        }
        val
      }
      gprime <- function(y, t){
        n <- length(y)
        val <- numeric(n)
        for(i in 1:n){
          val[i] <- ifelse(y[i] > 0, if(t[i] < 5) {return((exp(t[i]) * (-exp(exp(t[i])) + exp(t[i] + exp(t[i])) + 1)) / (exp(exp(t[i])) - 1)^2)}else{return(0)}, exp(t[i])) 
        }
        val
      }
    }
  }
  return(list(g = g, gprime = gprime))
}
