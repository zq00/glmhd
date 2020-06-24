#' Hinge function
#'
#' \eqn{\mathrm{hinge}(t) = \max(t,0)}
#' @param t Numeric. Input variable
#' @examples
#' hinge(1)
#' hinge(-1)
hinge <- function(t) max(t, 0)

