#' Methods for Linear Networks
#'
#' The function as.gn converts an object of class linnet to an object of
#' class gn
#'
#' @param object an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @return an object of class gn
#' @export
as.gn <- function(object){
  if (!inherits(object, "linnet")){
    object <- spatstat::as.linnet(object)
    if (!inherits(object, "linnet")){
      stop("Object must be of class 'linnet' or must be converable to an object of class 'linnet'")
    }
  }
  net <- list()
  net$vertices <- dplyr::tibble(x = object$vertices$x, y = object$vertices$y)

  object$d <- diag(object$dpath[object$from, object$to])
  A <- object$m*1
  V.deg <- as.numeric(table(c(object$from, object$to)))
  # vertex indices with degree 2
  ind.deg.2 <- which(V.deg == 2)

  # initialize
  ind.m.delete <- NULL
  i <- 0
  P <- list()

  # remove vertices with degree until none is left
  while(length(ind.deg.2) > 0){
    i <- i+1

    # find the two adjacent vertices to vertex with degree 2
    adj <- which(A[ind.deg.2[1], ] == 1)
    P[[i]] <- dplyr::tibble(from = c(adj[1], ind.deg.2[1]),
                     to = c(ind.deg.2[1], adj[2]),
                     m = NA,
                     length = NA)

    # go into the direction of the first adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[1], ]) == 2){
      adj.new <- which(A[adj[1], ] == 1)
      l <- adj.new[which(!is.element(adj.new, P[[i]]$from))]
      P[[i]] <- rbind(c(l, adj[1], NA, NA), P[[i]])
      adj[1] <- l
    }

    # go into the direction of the second adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[2], ]) == 2){
      adj.new <- which(A[adj[2], ] == 1)
      r <- adj.new[which(!is.element(adj.new, P[[i]]$to))]
      P[[i]] <- rbind(P[[i]], c(adj[2], r, NA, NA))
      adj[2] <- r
    }

    # save the line indices and their corresponding lengths which are removed from the network
    for (k in 1:nrow(P[[i]])) {
      P[[i]]$m[k] <- which(object$from == P[[i]]$from[k] & object$to == P[[i]]$to[k] | object$from == P[[i]]$to[k] & object$to == P[[i]]$from[k])
      P[[i]]$length[k] <- object$d[P[[i]]$m[k]]
    }

    # add line segments to the delete vector
    ind.m.delete <- unique(c(ind.m.delete, P[[i]]$m))

    # remove vertices from the current vector of vertices with degree 2
    ind.deg.2 <- setdiff(ind.deg.2, P[[i]]$to[1:(nrow(P[[i]])-1)])
    P[[i]]$m <- i
  }

  M <- object$lines$n - length(ind.m.delete) + length(P)
  curves <- vector("list", M)
  curves[1:length(P)] <- P
  k <- length(P)
  for (m in setdiff(1:object$lines$n, ind.m.delete)) {
    k <- k + 1
    curves[[k]] <- dplyr::tibble(from = object$from[m], to = object$to[m], m = k, length = object$d[m])
  }
  net$curves <- dplyr::bind_rows(curves)

  class(net) <- "gn"
  net
}

#' Methods for Geometric Networks
#'
#' The function as.gns converts an object of class gn to an object of
#' class gns
#'
#' @param object an object of class gn or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gn
#' @export
as.gns <- function(x, delta = NULL, h = NULL, r = 1){
  if (!inherits(x, c("gn", "linnet", "lpp"))){
    stop(paste("Object ", x, " can not be converted to an object of class gns"))
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
  e_to_v <- e_from_v <- vector("list", x$W)
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

  # incident lines to vertices
  ends_edges <- ends(x)
  for (v in 1:x$W) {
    vv <- which(x$A_v[v, ] == 1)
    if (length(vv) > 0){
      for (w in vv) {
        e_to_v[[v]] <- c(e_to_v[[v]], which(ends_edges$to == v & ends_edges$from == w))
        e_from_v[[v]] <- c(e_from_v[[v]], which(ends_edges$from == v & ends_edges$to == w))
      }
      e_to_v[[v]] <- sort(unique(e_to_v[[v]]))
      e_from_v[[v]] <- sort(unique(e_from_v[[v]]))
    }
  }
  x$splines <- list(delta = delta,
                    tau = tau,
                    J = J)
  x$bins <- list(h = h,
                 b = b,
                 z = z,
                 N = N)
  x$v_adj_e <- list(e_to_v = e_to_v,
                         e_from_v = e_from_v)
  B <- getB(x)
  x$B <- B
  x
}
