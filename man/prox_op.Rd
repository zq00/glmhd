% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prox_op.R
\name{prox_op}
\alias{prox_op}
\title{One-dimensional proximal operator}
\usage{
prox_op(f_prime, lambda, x)
}
\arguments{
\item{f_prime}{Function. Derivative of the function f.}

\item{lambda}{Numeric. Inverse penalty.}

\item{x}{Numeric. Input variable.}
}
\value{
A numeric value
    \deqn{\mathrm{prox}_{\lambda f}(x) = \mathrm{argmin}_x f(z) + \frac{1}{2\lambda}(z-x)^2.}
    It solves the equation
    \deqn{f'(z) + \frac{1}{\lambda}(z-x) = 0.}
}
\description{
One-dimensional proximal operator
}
\examples{
f <- function(x) 1/(1+exp(-x))
prox_op(f, 1, 1)
}
