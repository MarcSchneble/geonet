#' Confidence Bands of Smooth Terms
#'
#' \code{confidence_band} computes the lower and upper limits of smooth
#' terms fitted with \code{intensity_pspline}.
#'
#' @param theta The estimated coefficients which corresponds to the smooth
#' term.
#' @param V The covariance matrix of the estimated coefficients \code{theta}.
#' @param X The design matrix of the model which corresponds to the smooth
#' term.
#' @param q The quantile. Default to \code{q = 0.05} which corresponds to
#' 95% confidence bands.
#' @param R The number of replications in the simulation process.
#' @return A list of two vectors which contain the lower and the upper limits
#' of the confidence band.
#' @importFrom stats quantile
#' @importFrom mgcv rmvn

confidence_band <- function(theta, V, X, q = 0.05, R = 1000){
  mu_sim <- matrix(0, R, nrow(X))
  for (i in 1:R) {
    theta_sim <- rmvn(1, theta, V)
    mu_sim[i, ] <- as.vector(X%*%theta_sim)
  }
  lower <- upper <- rep(0, ncol(mu_sim))
  for (j in 1:ncol(mu_sim)) {
    lower[j] <- quantile(mu_sim[, j], probs = q/2)
    upper[j] <- quantile(mu_sim[, j], probs = 1 - q/2)
  }
  list(lower = lower, upper = upper)
}
