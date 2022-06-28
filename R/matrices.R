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
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}

incidence <- function(vertices, lins){
  e <- NULL
  W <- max(vertices$v, na.rm = TRUE)
  M <- max(lins$e, na.rm = TRUE)
  A <- Matrix::Matrix(0, W, M, sparse = TRUE)
  for (m in 1:M) {
    lins_m <- filter(lins, e == m)
    A[vertices$v[match(lins_m$a1[1], vertices$a)], m] <- -1
    A[vertices$v[match(utils::tail(lins_m$a2, 1), vertices$a)], m] <- 1
  }
  A
}

#' Design Matrix for Linear B-Splines on a Geometric Network
#'
#' \code{bspline_design} constructs the design matrix which represents
#' the (log-)baseline intensity on the geometric network.
#'
#' @param G An object of class \code{gn}.
#' @param knots A list which contains the knots on which the
#' B-splines are defined
#' @param bins list of which contains the mid points of the bins.
#' @return A sparse matrix design matrix of dimension N x J.
#' @importFrom splines splineDesign
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}

bspline_design <- function(G, knots, bins){
  # returns design matrix B with dimension N x J

  # design matrix for line segments
  B <- Matrix::Matrix(0, sum(bins$N), sum(knots$J) + G$W, sparse = TRUE)

  # line specific B-splines
  for (m in 1:G$M) {
    B[((cumsum(bins$N) - bins$N)[m] + 1):cumsum(bins$N)[m],
      ((cumsum(knots$J) - knots$J)[m] + 1):cumsum(knots$J)[m]] <-
      splineDesign(knots = knots$tau[[m]],
                            x = bins$z[[m]],
                            ord = 2,
                            outer.ok = TRUE,
                            sparse = TRUE)
  }

  # vertex specific B-splines
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      B[((cumsum(bins$N) - bins$N)[m] + 1):
          ((cumsum(bins$N) - bins$N)[m] +
             length(which(1 - bins$z[[m]]/knots$delta[m] > 0))), sum(knots$J) + v] <-
        (1 - bins$z[[m]]/knots$delta[m])[which(1 - bins$z[[m]]/knots$delta[m] > 0)]
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      B[(cumsum(bins$N)[m] - length(which(1 - (G$d[m] - bins$z[[m]])/knots$delta[m] > 0)) + 1):
          cumsum(bins$N)[m], sum(knots$J) + v] <-
        (1 - (G$d[m] - bins$z[[m]])/knots$delta[m])[which(1 - (G$d[m] - bins$z[[m]])/knots$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(bins$N)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}

#' B-Spline Design Matrix for Plotting
#'
#' @param X A point pattern on a geometric network (object of class \code{gnpp}).
#' @param df A data frame with points at which the fitted intensity should be
#' plotted.
#' @return A sparse design matrix.
#' @import dplyr
#' @importFrom splines splineDesign
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}

bspline_design_plot <- function(X, df){
  # returns design matrix for new data
  e <- NULL

  knots <- X$knots
  G <- X$network

  # design matrix for line segments
  B <- Matrix(0, nrow(df), sum(knots$J) + G$W, sparse = TRUE)
  N <- as.numeric(table(df$e))
  z <- vector("list", G$M)

  # line specific B-splines
  for (m in 1:G$M) {
    z[[m]] <- filter(df, e == m) %>% pull(z)
    B[((cumsum(N) - N)[m] + 1):cumsum(N)[m],
      ((cumsum(knots$J) - knots$J)[m] + 1):cumsum(knots$J)[m]] <-
      splineDesign(knots = knots$tau[[m]],
                   x = z[[m]], ord = 2, outer.ok = TRUE, sparse = TRUE)
  }


  # vertex specific B-splines
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      B[((cumsum(N) - N)[m] + 1):
          ((cumsum(N) - N)[m] +
             length(which(1 - (z[[m]])/knots$delta[m] > 0))), sum(knots$J) + v] <-
        (1 - z[[m]]/knots$delta[m])[which(1 - z[[m]]/knots$delta[m] > 0)]
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      B[(cumsum(N)[m] - length(which(1 - (G$d[m] - z[[m]])/knots$delta[m] > 0)) + 1):
          cumsum(N)[m], sum(knots$J) + v] <-
        (1 - (G$d[m] - z[[m]])/knots$delta[m])[which(1 - (G$d[m] - z[[m]])/knots$delta[m] > 0)]
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
#' \code{network_penalty} constructs the penalty matrix which relates to the
#' B-Splines created by \code{\link[geonet]{bspline_design}}.
#'
#' @param G A geometric network (object of class \code{gnpp}).
#' @param knots A list which contains the knots on which the
#' B-splines are defined
#' @param r The order of the penalty, default to first-order penalty (\code{r = 1}.
#' @return A sparse and square penalty matrix.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}

network_penalty <- function(G, knots, r){
  # returns the first or second penalty matrix K
  if (!r %in% c(1, 2)) stop("r must be either 1 or 2")
  # adjacency matrix of knots
  A_tau <- Matrix::Matrix(0, sum(knots$J) + G$W, sum(knots$J) + G$W, sparse = TRUE)

  # on each line
  for (m in 1:G$M) {
    if (knots$J[m] > 1){
      A_tau[((cumsum(knots$J) - knots$J)[m] + 2):cumsum(knots$J)[m],
            ((cumsum(knots$J) - knots$J)[m] + 1):(cumsum(knots$J)[m] - 1)] <-
        diag(knots$J[m] - 1)
    }
  }

  # around each vertex
  for (v in 1:G$W) {
    # left line ends
    for (m in which(G$incidence[v, ] == -1)) {
      A_tau[sum(knots$J) + v,
            (cumsum(knots$J) - knots$J)[m] + 1] <- 1
    }
    # right line ends
    for (m in which(G$incidence[v, ] == 1)) {
      A_tau[sum(knots$J) + v, cumsum(knots$J)[m]] <- 1
    }
  }

  # filling the upper diagonal
  A_tau <- A_tau + Matrix::t(A_tau)

  if (r == 1){
    # find for every spline function the adjacent spline functions
    adj <- Matrix::which(A_tau*lower.tri(A_tau) == 1, arr.ind = T)

    # initializing first order difference matrix
    D <- Matrix::Matrix(0, nrow(adj), sum(knots$J) + G$W, sparse = TRUE)

    for(i in 1:nrow(adj)){
      D[i, adj[i, 1]] <- 1
      D[i, adj[i, 2]] <- -1
    }
    return(Matrix::t(D)%*%D)
  }
  if (r == 2){
    # find for every spline function the adjacent spline functions
    S_A <- igraph::shortest.paths(igraph::graph_from_adjacency_matrix(A_tau))
    adj_1 <- which(S_A == 1, arr.ind = T)
    adj_2 <- which(S_A*lower.tri(S_A) == 2, arr.ind = T)

    # initializing second order difference matrix
    D <- Matrix::Matrix(0, nrow(adj_2), sum(knots$J) + G$W, sparse = TRUE)

    for (i in 1:nrow(adj_2)) {
      D[i, adj_2[i, 1]] <- D[i, adj_2[i, 2]] <- 1
      D[i, intersect(adj_1[which(adj_1[, 1] == adj_2[i, 1]), 2], adj_1[which(adj_1[, 1] == adj_2[i, 2]), 2])] <- -2
    }
    D <- Matrix::Matrix(D, sparse = TRUE)

    return(Matrix::t(D)%*%D)
  }
}
