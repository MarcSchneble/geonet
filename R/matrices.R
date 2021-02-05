#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param x an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
#' @import dplyr
#' @export

getIncidence <- function(x){
  e <- NULL
  incidence <- matrix(0, x$W, x$M)
  for (m in 1:x$M) {
    dat <- filter(x$lins, e == m)
    incidence[x$vertices$v[match(dat$v1[1], x$vertices$id)], m] <- -1
    incidence[x$vertices$v[match(utils::tail(dat$v2, 1), x$vertices$id)], m] <- 1
  }
  incidence
}

#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param X an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param P an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
#' @importFrom splines splineDesign
#' @export

getB <- function(X, P){
  # returns design matrix B with dimension N x J

  # design matrix for line segments
  B <- Matrix::Matrix(matrix(0, sum(P$bins$N), sum(P$splines$J) + X$W), sparse = TRUE)

  # line specific B-splines
  for (m in 1:X$M) {
    B[((cumsum(P$bins$N) - P$bins$N)[m] + 1):cumsum(P$bins$N)[m],
      ((cumsum(P$splines$J) - P$splines$J)[m] + 1):cumsum(P$splines$J)[m]] <-
      splineDesign(knots = P$splines$tau[[m]],
                            x = P$bins$z[[m]],
                            ord = 2,
                            outer.ok = TRUE,
                            sparse = TRUE)
  }

  # vertex specific B-splines
  for (v in 1:X$W) {
    # left line ends
    for (m in which(X$incidence[v, ] == -1)) {
      B[((cumsum(P$bins$N) - P$bins$N)[m] + 1):
          ((cumsum(P$bins$N) - P$bins$N)[m] +
             length(which(1 - P$bins$z[[m]]/P$splines$delta[m] > 0))), sum(P$splines$J) + v] <-
        (1 - P$bins$z[[m]]/P$splines$delta[m])[which(1 - P$bins$z[[m]]/P$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(X$incidence[v, ] == 1)) {
      B[(cumsum(P$bins$N)[m] - length(which(1 - (X$d[m] - P$bins$z[[m]])/P$splines$delta[m] > 0)) + 1):
          cumsum(P$bins$N)[m], sum(P$splines$J) + v] <-
        (1 - (X$d[m] - P$bins$z[[m]])/P$splines$delta[m])[which(1 - (X$d[m] - P$bins$z[[m]])/P$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(P$bins$N)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}

#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param X an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param df an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
#' @import dplyr
#' @importFrom splines splineDesign
#' @export

getBplot <- function(X, df){
  # returns design matrix for new data
  e <- NULL

  # design matrix for line segments
  B <- matrix(0, nrow(df), sum(X$model$P$splines$J) + X$W)
  N <- as.numeric(table(df$e))
  z <- vector("list", X$M)

  # line specific B-splines
  for (m in 1:X$M) {
    z[[m]] <- filter(df, e == m) %>% pull(z)
    B[((cumsum(N) - N)[m] + 1):cumsum(N)[m],
      ((cumsum(X$model$P$splines$J) - X$model$P$splines$J)[m] + 1):cumsum(X$model$P$splines$J)[m]] <-
      splineDesign(knots = X$model$P$splines$tau[[m]],
                            x = z[[m]], ord = 2, outer.ok = TRUE)
  }


  # vertex specific B-splines
  for (v in 1:X$W) {
    # left line ends
    for (m in which(X$incidence[v, ] == -1)) {
      B[((cumsum(N) - N)[m] + 1):
          ((cumsum(N) - N)[m] +
             length(which(1 - (z[[m]])/X$model$P$splines$delta[m] > 0))), sum(X$model$P$splines$J) + v] <-
        (1 - z[[m]]/X$model$P$splines$delta[m])[which(1 - z[[m]]/X$model$P$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(X$incidence[v, ] == 1)) {
      B[(cumsum(N)[m] - length(which(1 - (X$d[m] - z[[m]])/X$model$P$splines$delta[m] > 0)) + 1):
          cumsum(N)[m], sum(X$model$P$splines$J) + v] <-
        (1 - (X$d[m] - z[[m]])/X$model$P$splines$delta[m])[which(1 - (X$d[m] - z[[m]])/X$model$P$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (!all.equal(rowSums(B), rep(1, nrow(df)))){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}


#' Get the incidence matrix
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param X an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param P an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param r an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gnds
#' @export

getK <- function(X, P, r){
  # returns the first or second penalty matrix K

  # adjacency matrix of knots
  A_tau <- matrix(0, sum(P$splines$J) + X$W, sum(P$splines$J) + X$W)

  # on each line
  for (m in 1:X$M) {
    if (P$splines$J[m] > 1){
      A_tau[((cumsum(P$splines$J) - P$splines$J)[m] + 2):cumsum(P$splines$J)[m],
            ((cumsum(P$splines$J) - P$splines$J)[m] + 1):(cumsum(P$splines$J)[m] - 1)] <-
        diag(P$splines$J[m] - 1)
    }
  }

  # around each vertex
  for (v in 1:X$W) {
    # left line ends
    for (m in which(X$incidence[v, ] == -1)) {
      A_tau[sum(P$splines$J) + v,
            (cumsum(P$splines$J) - P$splines$J)[m] + 1] <- 1
    }
    # right line ends
    for (m in which(X$incidence[v, ] == 1)) {
      A_tau[sum(P$splines$J) + v, cumsum(P$splines$J)[m]] <- 1
    }
  }

  # filling the upper diagonal
  A_tau <- A_tau + t(A_tau)

  if (r == 1){
    # find for every spline function the adjacent spline functions
    adj <- which(A_tau*lower.tri(A_tau) == 1, arr.ind = T)

    # initalizing first order difference matrix
    D <- matrix(0, nrow(adj_2), sum(P$splines$J) + X$W)

    for(i in 1:nrow(adj)){
      D[i, adj[i, 1]] <- 1
      D[i, adj[i, 2]] <- -1
    }
    Matrix::t(D)%*%D
  }
  if (r == 2){
    # find for every spline function the adjacent spline functions
    S_A <- igraph::shortest.paths(igraph::graph_from_adjacency_matrix(A_tau))
    adj_1 <- which(S_A == 1, arr.ind = T)
    adj_2 <- which(S_A*lower.tri(S_A) == 2, arr.ind = T)

    # initalizing second order difference matrix
    D <- matrix(0, nrow(adj_2), sum(P$splines$J) + X$W)

    for (i in 1:nrow(adj_2)) {
      D[i, adj_2[i, 1]] <- D[i, adj_2[i, 2]] <- 1
      D[i, intersect(adj_1[which(adj_1[, 1] == adj_2[i, 1]), 2], adj_1[which(adj_1[, 1] == adj_2[i, 2]), 2])] <- -2
    }
    D <- Matrix::Matrix(D, sparse = TRUE)

    Matrix::t(D)%*%D
  }
}
