% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separable_proportion.R
\name{separable_proportion}
\alias{separable_proportion}
\title{Porportion of linearly separable subsamples}
\usage{
separable_proportion(X, Y, nn, B = 10)
}
\arguments{
\item{X}{Covariate matrix. Each row in \code{X} is one observation.}

\item{Y}{Response vector of \eqn{+1} and \eqn{-1} representing
the two classes. \code{Y} has the same length as the number of rows
in \code{X}.}

\item{nn}{Number of observations in each subsample.}

\item{B}{Numeric. Number of subsamples.}
}
\value{
Numeric. The proportion of separable subsamples.
}
\description{
This function randomly generate \code{B} subsamples of size \code{nn} and
returns the proportion of times that a subsample is linearly separable.
}
\examples{
\dontrun{
n <- 1000; p <- 400
X <- matrix(rnorm(n*p, 0, 1), n, p)
Y <- 2 * rbinom(n, 1, 0.5) - 1
separable_proportion(X, Y, nn = 600, B = 10)
}
}
\concept{is_separable}
