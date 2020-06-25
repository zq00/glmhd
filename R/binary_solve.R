#' Find solution of a 1d increasing function
#'
#' Binary search for zeros of a 1-d increasing function.
#' Values at the left and right side of the interval should be negative and positive respectively.
#'
#' @param f Function to minimize.
#' @param interval A vector of length 2. Interval that contains the zero.
#' @param eps Terminate when the size of interval is less than eps.
#' @param verbose Print progress if \code{TRUE}
#' @examples
#' \dontrun{
#' f <- function(t) t
#' binary_solve(f, c(-2, 2), eps = 0.001)
#' }
binary_solve <- function(f, interval, eps = 0.001, verbose = FALSE){
  l1 <- interval[1]; f1 <- f(l1)
  l2 <- interval[2]; f2 <- f(l2)
  l3 <- (l1 + l2) / 2; f3 <- f(l3)

  while(abs(l2 - l1) > eps){
    if(f3 > 0){
      l2 <- l3; f2 <- f3
    }else{
      l1 <- l3; f1 <- f3
    }
    l3 <- (l1 + l2) / 2; f3 <- f(l3)

    if(verbose)
      cat("Current interval is: [", l1, ",", l2, "]; Function values at endpoints are: [",  f1, ",", f2, "]\n")
  }
  return((l1 + l2) / 2)
}
