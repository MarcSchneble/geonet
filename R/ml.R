#' Maximum-Likelihood Estimation
#'
#' \code{logL}, \code{score} and \code{fisher} compute the log-likelihood,
#' the score-function and the Fisher-Information of the Poisson model.
#'
#' @param theta A vector of model coefficients.
#' @param design The design of the model.
#' @param rho A vector of smoothing parameters.
#' @return A scalar (\code{logL}), a vector of the same length as
#' theta (\code{score}) or a square matrix of the same dimension as
#' theta (\code{fisher}).
#' @export

score <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  score <- as.vector(Matrix::colSums((design$data$count - mu)*design$Z))
  for (a in 1:length(design$K)) {
    score[design$ind_smooths[[a]]] <-
      score[design$ind_smooths[[a]]] - rho[a]*design$K[[a]]%*%theta[design$ind_smooths[[a]]]
  }

  score
}

#' Maximum-Likelihood Estimation
#'
#' \code{logL}, \code{score} and \code{fisher} compute the log-likelihood,
#' the score-function and the Fisher-Information of the Poisson model.
#'
#' @param theta A vector of model coefficients.
#' @param design The design of the model.
#' @param rho A vector of smoothing parameters.
#' @return A scalar (\code{logL}), a vector of the same length as
#' theta (\code{score}) or a square matrix of the same dimension as
#' theta (\code{fisher}).
#' @export

fisher <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  Mu <- Matrix::Matrix(0, nrow(design$Z), nrow(design$Z))
  diag(Mu) <- mu
  fisher <- Matrix::t(design$Z)%*%Mu%*%design$Z
  for (a in 1:length(design$K)) {
    fisher[design$ind_smooths[[a]], design$ind_smooths[[a]]] <-
      fisher[design$ind_smooths[[a]], design$ind_smooths[[a]]] + rho[a]*design$K[[a]]
  }
  fisher
}

#' Maximum-Likelihood Estimation
#'
#' \code{logL}, \code{score} and \code{fisher} compute the log-likelihood,
#' the score-function and the Fisher-Information of the Poisson model.
#'
#' @param theta A vector of model coefficients.
#' @param design The design of the model.
#' @param rho A vector of smoothing parameters.
#' @param eps_theta asdasd
#' @return A scalar (\code{logL}), a vector of the same length as
#' theta (\code{score}) or a square matrix of the same dimension as
#' theta (\code{fisher}).
#' @export

scoring <- function(theta, design, rho, eps_theta = 1e-5){
  # perform iterative least squares estimation for a Poisson model with offset

  Delta_theta <- Inf
  it <- 0
  while (Delta_theta > eps_theta) {
    it <- it + 1
    theta_new <- as.vector(theta + Matrix::solve(fisher(theta, design, rho))%*%score(theta, design, rho))
    Delta_theta <- as.numeric(sqrt(t(theta - theta_new)%*%(theta - theta_new)))/
      as.numeric(sqrt(t(theta)%*%theta))
    theta <- as.vector(theta_new)
  }
  theta
}
