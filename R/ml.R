#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param theta an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param design an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param rho an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
#' @export

logL <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  logL <- sum(design$data$count*log(mu) - mu)
  for (a in 1:length(design$K)) {
    logL <- logL - 0.5*rho[a]*
      as.vector(theta[design$ind_smooths[[a]]]%*%design$K[[a]]%*%theta[design$ind_smooths[[a]]])
  }
  logL
}

#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param theta an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param design an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param rho an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
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

#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param theta an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param design an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param rho an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
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
