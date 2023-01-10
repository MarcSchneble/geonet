#' Maximum-Likelihood Estimation
#'
#' Scoring algorithm for maximum-likelihood estimation of a penalized Poisson
#' model while treating the smoothing parameters as fixed. Since the model
#' matrix \code{Z} when fitting a point process model on a geometric network is
#' very large with usually several millions of entries, \code{scoring} builds
#' an sparse representations of matrices in R.
#'
#' \code{scoring} performs the scoring algorithm for maximum-likelihood
#' estimation according to Fahrmeir et al. (2013). This algorithm is based
#' on the score-function and the Fisher-information of the log-likelihood.
#' \code{score} returns the score-function (the gradient of the log-likelihood)
#' and \code{fisher} returns the Fisher-information (negative Hessian of the
#' log-likelihood).
#'
#' @param theta An initial vector of model coefficients.
#' @param rho The current vector of smoothing parameters. For each smooth term,
#' including the baseline intensity of the network, one smoothing parameter
#' must be supplied.
#' @param data A data frame containing the data.
#' @param Z The (sparse) model matrix where the number of column must
#' correspond to the length of the vector of model coefficients \code{theta}.
#' @param K A (sparse) square penalty matrix of with the same dimension as
#' \code{theta}.
#' @param ind A list which contains the indices belonging to each smooth term
#' and the linear terms.
#' @param eps_theta The termination condition. If the relative change of the
#' norm of the model parameters is less than \code{eps_theta}, the scoring
#' algorithm terminates and returns the current vector of model parameters.
#' @return The maximum likelihood estimate for fixed smoothing parameters.
#' @references Fahrmeir, L., Kneib, T., Lang, S. and Marx, B. (2013).
#' Regression. Springer.
#' @keywords internal

scoring <- function(theta, rho, data, Z, K, ind, eps_theta = 1e-5){
  # perform iterative least squares estimation for a Poisson model with offset

  Delta_theta <- Inf
  it <- 0
  while (Delta_theta > eps_theta) {
    it <- it + 1
    theta_new <- as.vector(theta + Matrix::solve(fisher(theta, rho, data, Z, K, ind))%*%score(theta, rho, data, Z, K, ind))
    Delta_theta <- as.numeric(sqrt(t(theta - theta_new)%*%(theta - theta_new)))/
      as.numeric(sqrt(t(theta)%*%theta))
    theta <- as.vector(theta_new)
  }
  theta
}

#' @rdname scoring
#' @keywords internal

score <- function(theta, rho, data, Z, K, ind){
  mu <- exp(as.vector(Z%*%theta) + log(data$h) + log(data$offset))
  score <- as.vector(Matrix::colSums((data$count - mu)*Z)) -
    rep(c(rho, 0), lengths(ind))*K%*%theta
  score
}

#' @rdname scoring
#' @keywords internal

fisher <- function(theta, rho, data, Z, K, ind){
  mu <- exp(as.vector(Z%*%theta) + log(data$h) + log(data$offset))
  Mu <- Matrix::Matrix(0, nrow(Z), nrow(Z))
  diag(Mu) <- mu
  fisher <- Matrix::t(Z)%*%Mu%*%Z +
    rep(c(rho, 0), lengths(ind))*K
  fisher
}
