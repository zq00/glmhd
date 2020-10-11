#' Adjust the MLE in a binary regression
#'
#' This function computes the asymptotic bias and variance of the MLE of a binary regression
#' when both \eqn{n} and \eqn{p} are large,
#' assuming that covariates \eqn{X} are multivariate Gaussian.
#'
#' This function estimates coefficients of a binary regression by adjusting the MLE
#' based on the theory when \eqn{n} and \eqn{p} grows to infinity, and the covariates are multivariate Gaussian.
#' In this setting, the MLE of one variable \eqn{\beta_j} is
#' \deqn{
#' \hat{\beta}^{\mathrm{MLE}}_j - \alpha_\star \beta_j \approx \mathcal{N}(0, \sigma_\star^2 / \tau_j^2).
#' }
#' where \eqn{\tau_j^2} is the conditional variance of \eqn{x_j} given all the other variables.
#' This implies \eqn{\hat{\beta}_j^\mathrm{Adj} = \hat{\beta}^{MLE}_j / \alpha_\star} is unbiased for the coefficient \eqn{\beta_j}.
#' \code{adjust_binary} computes \eqn{\hat{\beta}_j^\mathrm{Adj}} and returns it in the vector \code{coef_adj}.
#'
#' @section Notes:
#' \itemize{
#'   \item The variables should be centered to have zero mean.
#'   \item The input to this function should be an object from the class \code{glm},
#'   and it should include \code{family} specifying the link function. Currently, the algorithm
#'   can handle logit and probit link.
#'   \item if gamma<0.0001, set gamma = 0
#' }
#'
#' @param glm_output An object from the class \code{glm},
#'     for example the output from the \code{\link[stats]{glm}} function
#' @param verbose If \code{TRUE}, print progress at every step.
#' @param echo If \code{TRUE}, return the input \code{glm_output}.
#' @param ... Additional arguments.
#' @return A list with elements
#' \describe{
#' \item{glm_output}{If \code{echo = TRUE}, returns the input \code{glm_output}.}
#' \item{gamma_hat}{Estimated signal strength \eqn{\gamma_0}.}
#' \item{param}{Estimated paramters \eqn{(\alpha_\star, \lambda_\star, \sigma_\star, b_\star)}}
#' \item{intercept}{Does the model contain an intercept?}
#' \item{tau_hat}{Estimated conditional standard deviation}
#' \item{coef_adj}{Adjusted MLE: this is \eqn{\hat{\beta}^{\mathrm{MLE}} / \alpha_\star}.}
#' \item{std_adj}{Estimated standard error of \eqn{\hat{\beta}_j}: \eqn{\sigma_\star / \tau_j}.}
#' \item{coef_unadj}{Unadjusted MLE.}
#' \item{std_unadj}{Unadjusted standard error. This is output from the \code{\link[stats]{glm}}.}
#' }
#' @include is_separable.R probe_frontier.R separable_proportion.R logistic_model.R
#'     integrate2_normal.R signal_strength.R prox_op.R find_param.R
#'     process_link.R h_eq.R hinge.R equation_binary.R
#' @export
#' @references
#' \emph{The Asymptotic Distribution of the MLE in High-dimensional Logistic Models: Arbitrary Covariance}, Qian Zhao, Pragya Sur and Emmanuel J. Candes, arXiv:2001.09351
#' @examples
#' \dontrun{
#' # Problem size
#' n <- 1000L
#' p <- 300L
#' # Generate data matrix
#' X <- matrix(rnorm(n*p, 0, 1), n, p) / sqrt(p)
#' # Sample parameters
#' beta <- rep(c(0, 1), p / 2) * 2
#' # Generate response
#' Y <- rbinom(n, 1, 1 / (1 + exp(- X %*% beta)))
#' fit <- glm(Y ~ X + 0, family = binomial, x = TRUE, y = TRUE)
#' adjusted_fit <- adjust_binary(fit)
#' # Adjusted MLE
#' head(adjusted_fit$coef_unadj)
#' }
#' @export
adjust_binary <- function(glm_output, verbose = TRUE, echo = TRUE, ...){
  # extract model matrix
  X <- glm_output$x
  Y <- (2 * glm_output$y) - 1
  # problem dimension
  n <- nrow(X); p <- ncol(X)
  kappa <- p/n
  # link functions
  link_fun <- process_link(glm_output$family)
  if(verbose) cat("Link function: ", glm_output$family$link, "\n")
  # is there an intercept in the model?
  has_intercept <- ifelse(names(glm_output$coef[1]) == "(Intercept)", TRUE, FALSE)
  # if the data is linearly separable, stop and print message
  if(is_separable(X,Y)) {cat("Data is separable. MLE does not exist!"); return(0)}
  # problem dimension that data becomes separable
  kappa_hat <- probe_frontier(X, Y, verbose = verbose, ...)

  # estimate signal strength
  p0 <- mean((Y+1)/2)
  signal_strength <- signal_strength(rho_prime = link_fun$rho_prime,
                                     kappa_hat,
                                     has_intercept, p0, verbose)
  # calculate parameters based on (kappa, gamma_hat)
  gamma_hat <- ifelse(signal_strength$gamma_hat > 1e-4, signal_strength$gamma_hat, 0)
  param <- find_param(
    rho_prime = link_fun$rho_prime,
    f_prime1 = link_fun$f_prime1,
    f_prime0 = link_fun$f_prime0,
    kappa,
    gamma = gamma_hat,
    intercept = has_intercept,
    beta0 = ifelse(!has_intercept, 0, signal_strength$b_hat),
    verbose = FALSE
  )

  if(has_intercept){names(param) <- c("alpha_s", "lambda_s", "sigma_s", "b_s")
  }else{ names(param) <- c("alpha_s", "lambda_s", "sigma_s")
  }

  if(verbose) {cat("Found parameters:\n"); print(param)}
  if(has_intercept){
    # Conditional standard deviation
    tau_hat <- 1 / sqrt(diag(solve(t(X[ , -1]) %*% X[ , -1]))) / sqrt(n) / sqrt(1 - kappa)
    # The adjusted mle and standard error
    coef_adj <-  glm_output$coefficients[-1] / param["alpha_s"]
    std_adj <- param["sigma_s"] / sqrt(p) / tau_hat
    return(
      list(
        glm_output = ifelse(echo, glm_output, NA),
        param = param,
        gamma_hat = gamma_hat,
        intercept = signal_strength$beta_hat * sign(glm_output$coef[1]),
        tau_hat = tau_hat,
        coef_adj = coef_adj,
        std_adj = std_adj,
        coef_unadj = glm_output$coefficients[-1],
        std_unadj = summary(glm_output)$coef[-1, 2]
      )
    )
  }else{
    # conditional standard deviation
    tau_hat <- 1 / sqrt(diag(solve(t(X) %*% X))) / sqrt(n) / sqrt(1 - kappa)
    # the adjusted mle and standard error
    coef_adj <- glm_output$coefficients / param["alpha_s"]
    std_adj <- param["sigma_s"] / sqrt(p) / tau_hat
    return(
      list(
        glm_output = ifelse(echo, glm_output, NA),
        gamma_hat = signal_strength$gamma_hat,
        param = param,
        intercept = NULL,
        tau_hat = tau_hat,
        coef_adj = coef_adj,
        std_adj = std_adj,
        coef_unadj = glm_output$coefficients,
        std_unadj = summary(glm_output)$coef[ , 2]
      )
    )
  }
}


