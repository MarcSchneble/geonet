#' Incidence matrix of a geometric network
#'
#' \code{incidence} constructs the incidence matrix of a geometric network
#' from the vertices and the curve segments.
#'
#' @param vertices A data frame containing the vertices of the geometric
#' network.
#' @param lins A data frame containing the curve segments of the geometric
#' network.
#' @return The incidence matrix of dimension \eqn{W} by \eqn{M}.
#' @import dplyr
#' @export

incidence <- function(vertices, lins){
  e <- NULL
  W <- max(vertices$v, na.rm = TRUE)
  M <- max(lins$e, na.rm = TRUE)
  A <- matrix(0, W, M)
  for (m in 1:M) {
    lins_m <- filter(lins, e == m)
    A[vertices$v[match(lins_m$v1[1], vertices$id)], m] <- -1
    A[vertices$v[match(utils::tail(lins_m$v2, 1), vertices$id)], m] <- 1
  }
  A
}

#' Design Matrix of Linear B-Splines on a Geometric Network
#'
#' \code{bspline_design} constructs the design matrix which represents
#' the baseline intensity on the geometric network.
#'
#' @param G A geometric network (object of class \code{gn}).
#' @param P B-splines on a geometric network.
#' @return A sparse matrix.
#' @importFrom splines splineDesign
#' @export

bspline_design <- function(G, P){
  # returns design matrix B with dimension N x J

  # design matrix for line segments
  B <- Matrix::Matrix(matrix(0, sum(P$bins$N), sum(P$splines$J) + G$W), sparse = TRUE)

  # line specific B-splines
  for (m in 1:G$M) {
    B[((cumsum(P$bins$N) - P$bins$N)[m] + 1):cumsum(P$bins$N)[m],
      ((cumsum(P$splines$J) - P$splines$J)[m] + 1):cumsum(P$splines$J)[m]] <-
      splineDesign(knots = P$splines$tau[[m]],
                            x = P$bins$z[[m]],
                            ord = 2,
                            outer.ok = TRUE,
                            sparse = TRUE)
  }

  # vertex specific B-splines
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      B[((cumsum(P$bins$N) - P$bins$N)[m] + 1):
          ((cumsum(P$bins$N) - P$bins$N)[m] +
             length(which(1 - P$bins$z[[m]]/P$splines$delta[m] > 0))), sum(P$splines$J) + v] <-
        (1 - P$bins$z[[m]]/P$splines$delta[m])[which(1 - P$bins$z[[m]]/P$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      B[(cumsum(P$bins$N)[m] - length(which(1 - (G$d[m] - P$bins$z[[m]])/P$splines$delta[m] > 0)) + 1):
          cumsum(P$bins$N)[m], sum(P$splines$J) + v] <-
        (1 - (G$d[m] - P$bins$z[[m]])/P$splines$delta[m])[which(1 - (G$d[m] - P$bins$z[[m]])/P$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(P$bins$N)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}

#' B-Spline Design for Plotting
#'
#'
#'
#' @param X A point pattern on a geometric network (object of class \code{gnpp}).
#' @param df A data frame with points at which the fitted intensity should be
#' plotted.
#' @return A sparse design matrix.
#' @import dplyr
#' @importFrom splines splineDesign
#' @export

getBplot <- function(X, df){
  # returns design matrix for new data
  e <- NULL

  P <- X$P
  G <- X$network

  # design matrix for line segments
  B <- Matrix(0, nrow(df), sum(P$splines$J) + G$W, sparse = TRUE)
  N <- as.numeric(table(df$e))
  z <- vector("list", G$M)

  # line specific B-splines
  for (m in 1:G$M) {
    z[[m]] <- filter(df, e == m) %>% pull(z)
    B[((cumsum(N) - N)[m] + 1):cumsum(N)[m],
      ((cumsum(P$splines$J) - P$splines$J)[m] + 1):cumsum(P$splines$J)[m]] <-
      splineDesign(knots = P$splines$tau[[m]],
                   x = z[[m]], ord = 2, outer.ok = TRUE, sparse = TRUE)
  }


  # vertex specific B-splines
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      B[((cumsum(N) - N)[m] + 1):
          ((cumsum(N) - N)[m] +
             length(which(1 - (z[[m]])/P$splines$delta[m] > 0))), sum(P$splines$J) + v] <-
        (1 - z[[m]]/P$splines$delta[m])[which(1 - z[[m]]/P$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      B[(cumsum(N)[m] - length(which(1 - (G$d[m] - z[[m]])/P$splines$delta[m] > 0)) + 1):
          cumsum(N)[m], sum(P$splines$J) + v] <-
        (1 - (G$d[m] - z[[m]])/P$splines$delta[m])[which(1 - (G$d[m] - z[[m]])/P$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (!all.equal(Matrix::rowSums(B), rep(1, nrow(df)))){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}


#' Penalty Matrix of a Geometric Network
#'
#' \code{penalty_network} constructs the penalty matrix which relates to the B-Spline
#' design matrix created by \code{\link[geonet]{bspline_design}}.
#'
#' @param G A geometric network (object of class \code{gnpp}).
#' @param P B-splines on a geometric network.
#' @param r The order of the penalty, default to first-order penalty (\code{r = 1}.
#' @return A sparse and square penalty matrix.
#' @export

penalty_network <- function(G, P, r){
  # returns the first or second penalty matrix K
  if (!r %in% c(1, 2)) stop("r must be either 1 or 2")
  # adjacency matrix of knots
  A_tau <- matrix(0, sum(P$splines$J) + G$W, sum(P$splines$J) + G$W)

  # on each line
  for (m in 1:G$M) {
    if (P$splines$J[m] > 1){
      A_tau[((cumsum(P$splines$J) - P$splines$J)[m] + 2):cumsum(P$splines$J)[m],
            ((cumsum(P$splines$J) - P$splines$J)[m] + 1):(cumsum(P$splines$J)[m] - 1)] <-
        diag(P$splines$J[m] - 1)
    }
  }

  # around each vertex
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      A_tau[sum(P$splines$J) + v,
            (cumsum(P$splines$J) - P$splines$J)[m] + 1] <- 1
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      A_tau[sum(P$splines$J) + v, cumsum(P$splines$J)[m]] <- 1
    }
  }

  # filling the upper diagonal
  A_tau <- A_tau + t(A_tau)

  if (r == 1){
    # find for every spline function the adjacent spline functions
    adj <- which(A_tau*lower.tri(A_tau) == 1, arr.ind = T)

    # initalizing first order difference matrix
    D <- matrix(0, nrow(adj_2), sum(P$splines$J) + G$W)

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
    D <- matrix(0, nrow(adj_2), sum(P$splines$J) + G$W)

    for (i in 1:nrow(adj_2)) {
      D[i, adj_2[i, 1]] <- D[i, adj_2[i, 2]] <- 1
      D[i, intersect(adj_1[which(adj_1[, 1] == adj_2[i, 1]), 2], adj_1[which(adj_1[, 1] == adj_2[i, 2]), 2])] <- -2
    }
    D <- Matrix::Matrix(D, sparse = TRUE)

    Matrix::t(D)%*%D
  }
}
