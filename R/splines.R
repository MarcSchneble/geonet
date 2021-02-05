#' Methods for Geometric Networks
#'
#' The function as.gnds converts an object of class gnpp to an object of
#' class gnds
#'
#' @param X an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @return an object of class gnds
#' @export

Pspline <- function(X, delta = NULL, h = NULL, r = 1){
  if (!inherits(X, c("gn", "gnpp"))){
    stop(paste("Object ", deparse(quote(X)), " can not be converted to an object of class gns"))
  }
  if (inherits(X, c("linnet", "lpp"))) xX <- as.gn(X)
  if (is.null(delta)) delta <- min(X$d/2)
  if (is.null(h)) h <- delta/2
  if (!r %in% c(1, 2)) stop("r must be either 1 or 2")
  # line specific knot distances
  delta <- X$d*(delta > X$d) + delta*(delta <= X$d)
  delta <- pmin(X$d/floor(X$d/delta)*(X$d/delta - floor(X$d/delta) < 0.5) +
                  X$d/ceiling(X$d/delta)*(X$d/delta - floor(X$d/delta) >= 0.5), X$d/2)
  # ensure that h <= delta
  h <- min(h, min(delta))
  # line specific bin widths
  h <- X$d/floor(X$d/h)*(X$d/h - floor(X$d/h) < 0.5) +
    X$d/ceiling(X$d/h)*(X$d/h - floor(X$d/h) >= 0.5)
  # initializing...
  tau <- b <- z <- vector("list", X$M)
  N <- J <- rep(0, X$M)

  # do for every line segment
  for (m in 1:X$M) {

    # knot sequences tau
    tau[[m]] <- seq(0, X$d[m], delta[m])

    # bin boundaries b
    b[[m]] <- seq(0, X$d[m], h[m])

    # characterization of bins by midpoints z
    z[[m]] <- (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2

    # total count of bins in the geometric network
    N[m] <- length(z[[m]])

    # count of linear B-splines on line segment
    J[m] <- length(tau[[m]]) - 2
  }
  P <- list()
  P$splines <- list(delta = delta,
                    tau = tau,
                    J = J)
  P$bins <- list(h = h,
                 b = b,
                 z = z,
                 N = N)
  P$B <- getB(X, P)
  P$K <- getK(X, P, r)
  P
}
