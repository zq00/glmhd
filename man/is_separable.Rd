% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_separable.R
\name{is_separable}
\alias{is_separable}
\title{Are the two classes linearly separable?}
\usage{
is_separable(X, Y, add_intercept = FALSE)
}
\arguments{
\item{X}{Covariate matrix. Each row of \code{X} corresponds to one observation.}

\item{Y}{A vector of \eqn{+1} and \eqn{-1} representing the two classes.
\code{Y} has the same length as the number of rows in \code{X}.}

\item{add_intercept}{If \code{TRUE}, add an intercept to the matrix \code{X}.
Set to \code{FALSE} by default.}
}
\value{
Returns 1 if the data is separable and 0 otherwise.
}
\description{
The function \code{is_separable} determines if the two classes \eqn{Y=1}
and \eqn{Y=-1} are linearly separable using columns of \eqn{X}.
}
\details{
The two classes are separable if there exists a non-zero vector \eqn{b}
that satisfies \deqn{y_i x_i^\top b \leq 0} for every observation \eqn{i}.
}
\examples{
\dontrun{
n <- 1000; p <- 400
X <- matrix(rnorm(n*p, 0, 1), n, p)
Y <- 2 * rbinom(n, 1, 0.5) - 1
is_separable(X, Y, add_intercept = TRUE)
}
}
\references{
\emph{The phase transition for the existence of the maximum likelihood estimate in high-dimensional logistic regression}
Emmanuel J. Candes and Pragya Sur, Ann. Statist., Volume 48, Number 1 (2020), 27-42.

\emph{Linear programming algorithms for detecting separated data in binary logistic regression models}, Kjell Konis, Ph.D. thesis, Univ. Oxford.
}
\seealso{
\link[safeBinaryRegression]{glm}
}
\keyword{internal}
