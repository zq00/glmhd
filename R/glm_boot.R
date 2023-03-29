#' Resized bootstrap method for a GLM
#' 
#' Estimates the MLE distribution of a GLM using the resized bootstrap method
#' 
#' @param glm_fit A glm object returned by [stats:glm] function. It should contain the covariates x and responses y, i.e., the output of glm(.., x = T, y = T).
#' @param simulate_fun A function to simulate new responses (the inputs are an observation matrix and a vector of coefficients).
#' @param s_interval A numeric value of increment of the sequence of shrinkage factors \eqn{s}. If \eqn{s_interval = 0.2}, then the 
#'   list of shrinkage factors would be \eqn{(0, 0.2, 0.4, ..., 1)}. 
#' @param b_var A numeric value of the number of parametric bootstrap samples at each s to estimate signal strength parameter \eqn{\gamma}.
#' @param b_boot A numeric value of the number of bootstrap samples to estimate the bias and variance of the MLE.
#' @param robust_est If \code{TRUE}, use robust estimator of the bias and std.dev.
#' @param filename filename If a file name is provided, then save the plot of \eqn{\hat{\eta}(s)} versus shrinkage factors \eqn{\gamma} to \code{filename} (see also [estimate_gamma]).
#' @param verbose Print progress if \code{TRUE}.
#' @return 
#' \describe{
#' \item{glm_fit}{The input GLM object.}
#' \item{beta_s}{\eqn{\beta_s = s_\star\cdot \hat{\beta}} satisfies \eqn{\var(X^\top \beta_s)\approx \gamma^2}.}
#' \item{gamma_hat}{The estimated signal strength parameter \eqn{\gamma}.}
#' \item{alpha}{A numeric value of the estimated inflation of the MLE, i.e., \eqn{\hat{\beta}/\alpha} is approximately unbiased of the true model coef. }
#' \item{sd}{A numeric value of the estimated std.dev. of the MLE.}
#' \item{boot_sample}{A matrix of size {p}*\code{b_boot} (p is the number of variables) of the bootstrap MLE. }
#' }
#' @importFrom robustbase Qn 
#' @examples
#' Problem size
#' n <- 1000L
#' p <- 300L
#' # Generate data matrix
#' nu <-  10
#' chi <- rchisq(n, df = nu) / (nu - 2)
#' X <- matrix(rnorm(n * p, 0, 1), n, p) / sqrt(chi) / sqrt(n)
#' # Sample parameters
#' beta <- rep(c(0, 1), p / 2) * 2
#' # Generate response
#' Y <- rbinom(n, 1, 1 / (1 + exp(- X %*% beta)))
#' fit <- glm(Y ~ X + 0, family = binomial, x = TRUE, y = TRUE)
#' adjusted_fit <- glm_boot(fit)
#' # Estimated inflation
#' adjusted_fit$alpha
#' # Esimtaed std.dev.
#' head(adjusted_fit$sd)
#' @export
#' @export
glm_boot <- function(glm_fit, simulate_fun = NULL, s_interval = 0.02, b_var = 5, b_boot = 100, robust_est = FALSE, verbose = TRUE, filename = NA){
  # 1. extract data and MLE
  family <- glm_fit$family # family and link 
  if(is.null(simulate_fun)){
    family$simulate_fun <- get_simulate_fun(family)
  }else{
    family$simulate_fun <- simulate_fun
  }
  
  X <- glm_fit$x; Y <- glm_fit$y # covariate matrix contains a first column of 1s if the model contains an intercept
  n <- nrow(X); p <- ncol(X) 
  beta_hat <- glm_fit$coef
  
  # 2. Estimate gamma and beta0
  eta_obs <- estimate_eta(X, Y, beta_hat, family) # estimates eta from the obs. 
  if(verbose) cat("Observed eta = ", eta_obs, "\n")
  if(eta_obs > 1000) stop("Error: estimated eta is too large!") 
    
  # Estimate eta_hat at a sequence of shrinkage factors
  s_seq <- seq(0, 1, by = s_interval); ns <- length(s_seq)
  eta_hat <- matrix(0, ns, b_var); i <- 0
  for(s in s_seq){ 
    new_val <- estimate_variance(X, beta_hat * s, family, b_var)
  
    if((is.numeric(new_val) && length(new_val) == 1 && new_val == -1) || mean(new_val) > 1.5 * eta_obs) break; # stop criteria
    i <- i+1; eta_hat[i, ] <- new_val
    if(verbose){if(i %% 2 == 0){cat(s_seq[i], "\t Estimated eta is ", mean(eta_hat[i, ]),"\n") }}
  }
  if(i == 1 || i == 0) {s <- 0; sol <- list(gamma_hat = 0, s_hat=0)}else{
    s_seq <- s_seq[1:i]; eta_hat <- eta_hat[1:i, ]
    # find solutions s_seq, eta_hat, eta_obs, sd_obs, verbose = T, filename = NULL
    sol <- estimate_gamma( s_seq, eta_hat, eta_obs, sd(X%*% beta_hat), verbose = verbose, filename = filename)
  }
  s_hat <- sol$s_hat
  beta_s <- beta_hat * s_hat; 
  if(verbose){cat("Estimated gamma is", sol$gamma_hat , "\n")}
  
  # 3. Using bootstrap to estimate the bias and variance
  mle_boot <- bootglm(X, beta_s, family, b_boot, verbose)
  
  # 4. estimate alpha and sigma
  if(robust_est){
    sd_boot <- apply(mle_boot, 1, function(t) Qn(t, constant = 2.2219, finite.corr = F))
  }else{
    sd_boot <- apply(mle_boot, 1, sd)
  }
  
  if(s_hat == 0){ 
    alpha_boot <- 1
    cat("No signal! \n")
  }else{ 
    if(robust_est){
      mean_ap <- apply(mle_boot, 1, median)
    }else{
      mean_ap <- rowMeans(mle_boot)
    }
    alpha_boot <- lm(mean_ap ~ beta_s + 0, weights = 1/sd_boot^2)$coef
  }
  
  return(list(glm_fit = glm_fit, beta_s = beta_s, gamma_hat = sol$gamma_hat, alpha = alpha_boot, sd = sd_boot, boot_sample = mle_boot))
}

