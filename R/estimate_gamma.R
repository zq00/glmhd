#' Estimate the signal strength 
#' 
#' Estimating gamma using the SLOE estimator and parametric bootstrap
#' 
#' We estimate \eqn{\gamma} by the standard deviation of \eqn{sd(x_i^\top \beta(s_{\star}))}, 
#'   where \eqn{\beta(s_\star) = s_\star \cdot \hat{\beta}}. We use the following relationship: if 
#'   \eqn{sd(x_i^\top \beta(s_{\star}))\approx \gamma}, then \eqn{\hat{\eta}(s_\star) \approx \hat{\eta}}.
#'   In this equation, \eqn{\hat{\eta}} is the estimated \eqn{\hat{\eta}} from the observations,
#'   \eqn{\eta(s)} is \eqn{\eta} when the model coefficient is \eqn{\beta(s) = s\cdot\hat{\beta}}. We 
#'   estimate \eqn{\eta(s)} by parametric bootstrap, fixing the covariates at the observed values and setting the model 
#'   coefficients as \eqn{\beta(s)} (see [estimate_variance] on how to estimate \eqn{\eta(s)}). 
#'   
#' We pick a sequence of shrinkage factors \eqn{s}, and then compute \eqn{\hat{\eta}(s)} for each of them.
#'   Then, we fit a LOESS curve of \eqn{\eta(s)} as a function of \eqn{s} using the sequence of \eqn{s} and 
#'   the estimated \eqn{\hat{\eta}(s)} at each bootstrap sample. Finally, we choose the shrinkage factor \eqn{s_\star} on the 
#'   curve such that the fitted value is equal to the observed \eqn{\hat{\eta}}. 
#'   
#' The estimated \eqn{\hat{\gamma}} is \eqn{sd(x_i^\top \beta(s_{\star}))}. 
#'  
#' @param s_seq A sequence of shrinkage factors s.
#' @param eta_hat A matrix. The number of rows is equal to the length of s_seq, the number of columns is equal to the number of parametric bootstrap samples at each s.
#' @param eta_obs \eqn{\hat{\eta}} computed using observed samples.
#' @param sd_obs Observed standard deviation of the linear predictors evaluated at \eqn{\hat{\beta}}, i.e., \eqn{sd(x_i^\top \hat{\beta})}.
#' @param verbose Plot \eqn{\gamma} versus shrinkage factors \eqn{s} if \code{TRUE}.
#' @param filename If a file name is provided, then save the plot of \eqn{\hat{\eta}(s)} versus \eqn{\gamma} to \code{filename}.

#' @return 
#' \describe{
#' \item{s_hat}{A numeric value of the estimated shrinkage factor \eqn{s} that satisfies
#' \eqn{\hat{s} \mathrm{sd}(X^\top \hat{\beta}) = \hat{\gamma}.}
#' }
#' \item{gamma_hat}{A numeric value of the estimated signal strength.}
#' }
#' @importFrom stats loess 
#' @importFrom stats predict 
#' @importFrom stats sd
#' @importFrom stats median 
#' @importFrom stats lm 
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics mtext
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @export
estimate_gamma <- function(s_seq, eta_hat, eta_obs, sd_obs, verbose = T, filename = NA){
  # fit a smooth loess curve 
  data <- data.frame(cbind(rep(s_seq, time = ncol(eta_hat)), as.vector(eta_hat)))
  colnames(data) <- c("s", "val")
  curve <- loess(val ~ s, data = data)
  s_new <- seq(min(s_seq), max(s_seq), by = 0.001)
  eta_new <- predict(curve, s_new) # estimated sd_hat on the smoothed loess curve
  diff <- abs(eta_new - eta_obs)
  s_sol <- s_new[which.min(diff)]
  gamma_hat <- sd_obs * s_sol # estimated gamma_hat
  
  if(verbose){ # plot eta(s) versus gamma(s)
    if(!is.na(filename)){png(filename, width = 500, height = 400)} # if input a file location, save plot to the file location
    gamma <- sd_obs * s_new
    plot(gamma, eta_new, type = "l", ylim = range(eta_new), 
         xlab = expression(gamma), ylab = expression(paste(hat(eta))))
    points( sd_obs * data$s, data$val, pch = 16, cex = 0.5)
    abline(h = eta_obs, lty = "dotted")
    mtext(text = paste0("gamma_hat = ", round(gamma_hat, 2)), line = 0.5)
    if(!is.na(filename)) dev.off()
  }
  return(list(s_hat = s_sol, gamma_hat = gamma_hat))
}
