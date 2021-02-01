logL <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  logL <- sum(design$data$count*log(mu) - mu)
  for (a in 1:length(design$K)) {
    logL <- logL - 0.5*rho[a]*
      as.vector(theta[design$ind_smooths[[a]]]%*%design$K[[a]]%*%theta[design$ind_smooths[[a]]])
  }
  return(logL)
}

score <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  score <- as.vector(Matrix::colSums((design$data$count - mu)*design$Z))
  for (a in 1:length(design$K)) {
    score[design$ind_smooths[[a]]] <-
      score[design$ind_smooths[[a]]] - rho[a]*design$K[[a]]%*%theta[design$ind_smooths[[a]]]
  }

  return(score)
}

fisher <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  Mu <- Matrix::Matrix(0, nrow(design$Z), nrow(design$Z))
  diag(Mu) <- mu
  fisher <- Matrix::t(design$Z)%*%Mu%*%design$Z
  for (a in 1:length(design$K)) {
    fisher[design$ind_smooths[[a]], design$ind_smooths[[a]]] <-
      fisher[design$ind_smooths[[a]], design$ind_smooths[[a]]] + rho[a]*design$K[[a]]
  }
  return(fisher)
}
