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
  if (!inherits(X, c("gnpp"))){
    stop(paste("Object ", deparse(quote(X)), " can not be converted to an object of class gns"))
  }
  if (inherits(X, c("linnet", "lpp"))) X <- as.gnpp(X)
  if (is.null(delta)) delta <- min(G$d)/2
  if (is.null(h)) h <- delta/2
  if (!r %in% c(1, 2)) stop("r must be either 1 or 2")
  # line specific knot distances
  G <- as.gn(X)
  delta <- G$d*(delta > G$d) + delta*(delta <= G$d)
  delta <- pmin(G$d/floor(G$d/delta)*(G$d/delta - floor(G$d/delta) < 0.5) +
                  G$d/ceiling(G$d/delta)*(G$d/delta - floor(G$d/delta) >= 0.5), G$d/2)
  # ensure that h <= delta
  h <- min(h, min(delta))
  # line specific bin widths
  h <- G$d/floor(G$d/h)*(G$d/h - floor(G$d/h) < 0.5) +
    G$d/ceiling(G$d/h)*(G$d/h - floor(G$d/h) >= 0.5)
  # initializing...
  tau <- b <- z <- vector("list", G$M)
  N <- J <- rep(0, G$M)

  # do for every line segment
  for (m in 1:G$M) {

    # knot sequences tau
    tau[[m]] <- seq(0, G$d[m], delta[m])

    # bin boundaries b
    b[[m]] <- seq(0, G$d[m], h[m])

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
  P$B <- getB(G, P)
  P$K <- getK(G, P, r)
  P
}
