#' Methods for Geometric Networks
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param x an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @return an object of class gnds
#' @export
addSpline <- function(x, delta = NULL, h = NULL, r = 1){
  if (!inherits(x, c("gn", "gnpp", "linnet", "lpp"))){
    stop(paste("Object ", deparse(quote(x)), " can not be converted to an object of class gns"))
  }
  if (inherits(x, c("linnet", "lpp"))) x <- as.gn(x)
  if (is.null(delta)) delta <- min(x$d/2)
  if (is.null(h)) h <- delta/2
  if (!r %in% c(1, 2)) stop("r must be either 1 or 2")
  # line specific knot distances
  delta <- x$d*(delta > x$d) + delta*(delta <= x$d)
  delta <- pmin(x$d/floor(x$d/delta)*(x$d/delta - floor(x$d/delta) < 0.5) +
                  x$d/ceiling(x$d/delta)*(x$d/delta - floor(x$d/delta) >= 0.5), x$d/2)
  # ensure that h <= delta
  h <- min(h, min(delta))
  # line specific bin widths
  h <- x$d/floor(x$d/h)*(x$d/h - floor(x$d/h) < 0.5) +
    x$d/ceiling(x$d/h)*(x$d/h - floor(x$d/h) >= 0.5)
  # initializing...
  tau <- b <- z <- vector("list", x$M)
  N <- J <- rep(0, x$M)

  # do for every line segment
  for (m in 1:x$M) {

    # knot sequences tau
    tau[[m]] <- seq(0, x$d[m], delta[m])

    # bin boundaries b
    b[[m]] <- seq(0, x$d[m], h[m])

    # characterization of bins by midpoints z
    z[[m]] <- (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2

    # total count of bins in the geometric network
    N[m] <- length(z[[m]])

    # count of linear B-splines on line segment
    J[m] <- length(tau[[m]]) - 2
  }
  x$splines <- list(delta = delta,
                    tau = tau,
                    J = J)
  x$bins <- list(h = h,
                 b = b,
                 z = z,
                 N = N)
  x$B <- getB(x)
  x$K <- getK(x, r)
  class(x) <- c(class(x), "gnpp")
  x
}

#' B-Spline Model Matrix B
#'
#' The function calculates the B-spline model matrix B
#'
#' @param x an object of class gn augmented by
#' @return model matrix B
#' @export

getB <- function(x){
  # returns design matrix B with dimension N x J

  # design matrix for line segments
  B <- Matrix::Matrix(matrix(0, sum(x$bins$N), sum(x$splines$J) + x$W), sparse = TRUE)

  # line specific B-splines
  for (m in 1:x$M) {
    B[((cumsum(x$bins$N) - x$bins$N)[m] + 1):cumsum(x$bins$N)[m],
      ((cumsum(x$splines$J) - x$splines$J)[m] + 1):cumsum(x$splines$J)[m]] <-
      splines::splineDesign(knots = x$splines$tau[[m]],
                            x = x$bins$z[[m]],
                            ord = 2,
                            outer.ok = TRUE,
                            sparse = TRUE)
  }

  # vertex specific B-splines
  for (v in 1:x$W) {
    # left line ends
    for (m in which(x$incidence[v, ] == -1)) {
      B[((cumsum(x$bins$N) - x$bins$N)[m] + 1):
          ((cumsum(x$bins$N) - x$bins$N)[m] +
             length(which(1 - x$bins$z[[m]]/x$splines$delta[m] > 0))), sum(x$splines$J) + v] <-
        (1 - x$bins$z[[m]]/x$splines$delta[m])[which(1 - x$bins$z[[m]]/x$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(x$incidence[v, ] == 1)) {
      B[(cumsum(x$bins$N)[m] - length(which(1 - (x$d[m] - x$bins$z[[m]])/x$splines$delta[m] > 0)) + 1):
          cumsum(x$bins$N)[m], sum(x$splines$J) + v] <-
        (1 - (x$d[m] - x$bins$z[[m]])/x$splines$delta[m])[which(1 - (x$d[m] - x$bins$z[[m]])/x$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(x$bins$N)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}

#' Penalty Matrix K
#'
#' The function calculates the penalty matrix L
#'
#' @param x an object of class gn augmented by
#' @param r the order of the penalty
#' @return penalty matrix K
#' @export

getK <- function(x, r){
  # returns the first or second penalty matrix K

  # adjacency matrix of knots
  A_tau <- matrix(0, sum(x$splines$J) + x$W, sum(x$splines$J) + x$W)

  # on each line
  for (m in 1:x$M) {
    if (x$splines$J[m] > 1){
      A_tau[((cumsum(x$splines$J) - x$splines$J)[m] + 2):cumsum(x$splines$J)[m],
            ((cumsum(x$splines$J) - x$splines$J)[m] + 1):(cumsum(x$splines$J)[m] - 1)] <-
        diag(x$splines$J[m] - 1)
    }
  }

  # around each vertex
  for (v in 1:x$W) {
    # left line ends
    for (m in which(x$incidence[v, ] == -1)) {
      A_tau[sum(x$splines$J) + v,
            (cumsum(x$splines$J) - x$splines$J)[m] + 1] <- 1
    }
    # right line ends
    for (m in which(x$incidence[v, ] == 1)) {
      A_tau[sum(x$splines$J) + v, cumsum(x$splines$J)[m]] <- 1
    }
  }

  # filling the upper diagonal
  A_tau <- A_tau + t(A_tau)

  if (r == 1){
    # find for every spline function the adjacent spline functions
    adj <- which(A_tau*lower.tri(A_tau) == 1, arr.ind = T)

    # initalizing first order difference matrix
    D <- Matrix::Matrix(matrix(0, nrow(adj), sum(x$splines$J) + x$W), sparse = TRUE)

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

    # initalizing second order difference matrix
    D <- Matrix::Matrix(matrix(0, nrow(adj_2), sum(x$splines$J) + x$W), sparse = TRUE)

    for (i in 1:nrow(adj_2)) {
      D[i, adj_2[i, 1]] <- D[i, adj_2[i, 2]] <- 1
      D[i, intersect(adj_1[which(adj_1[, 1] == adj_2[i, 1]), 2], adj_1[which(adj_1[, 1] == adj_2[i, 2]), 2])] <- -2
    }
    return(Matrix::t(D)%*%D)
  }
}

getBplot <- function(x, df){
  # returns design matrix for new data

  # design matrix for line segments
  B <- matrix(0, nrow(df), sum(x$splines$J) + x$W)
  N <- as.numeric(table(df$e))
  z <- vector("list", x$M)

  # line specific B-splines
  for (m in 1:x$M) {
    z[[m]] <- dplyr::filter(df, e == m) %>% dplyr::pull(z)
    B[((cumsum(N) - N)[m] + 1):cumsum(N)[m],
      ((cumsum(x$splines$J) - x$splines$J)[m] + 1):cumsum(x$splines$J)[m]] <-
      splines::splineDesign(knots = x$splines$tau[[m]],
                   x = z[[m]], ord = 2, outer.ok = TRUE)
  }


  # vertex specific B-splines
  for (v in 1:x$W) {
    # left line ends
    for (m in which(x$incidence[v, ] == -1)) {
      B[((cumsum(N) - N)[m] + 1):
          ((cumsum(N) - N)[m] +
             length(which(1 - (z[[m]])/x$splines$delta[m] > 0))), sum(x$splines$J) + v] <-
        (1 - z[[m]]/x$splines$delta[m])[which(1 - z[[m]]/x$splines$delta[m] > 0)]
    }
    # right line ends
    for (m in which(x$incidence[v, ] == 1)) {
      B[(cumsum(N)[m] - length(which(1 - (x$d[m] - z[[m]])/x$splines$delta[m] > 0)) + 1):
          cumsum(N)[m], sum(x$splines$J) + v] <-
        (1 - (x$d[m] - z[[m]])/x$splines$delta[m])[which(1 - (x$d[m] - z[[m]])/x$splines$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (!all.equal(rowSums(B), rep(1, nrow(df)))){
    stop("Error! Rowsums of B are not equal to one!")
  }
  B
}
