#' Methods for Linear Networks
#'
#' The function as.gn converts an object of class linnet to an object of
#' class gn
#'
#' @param object an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @return an object of class gn
#' @export
as.gn <- function(L){
  if (!inherits(L, "linnet")){
    L <- spatstat::as.linnet(L)
    if (!inherits(L, "linnet")){
      stop("Object must be of class 'linnet' or must be converable to an object of class 'linnet'")
    }
  }
  G <- list(
    vertices = dplyr::tibble(id = 1:L$vertices$n,
                             v = NA,
                             x = L$vertices$x,
                             y = L$vertices$y),
    lins = dplyr::tibble(id = 1:L$lines$n,
                         e = NA,
                         v1 = L$from,
                         v2 = L$to,
                         v1_x = L$lines$ends$x0,
                         v1_y = L$lines$ends$y0,
                         v2_x = L$lines$ends$x1,
                         v2_y = L$lines$ends$y1,
                         length = diag(L$dpath[L$from, L$to])),
    adjacency = NULL, incidence = NULL,
    d = NULL, q = 2, W = NULL, M = NULL
  )
  # adjacency matrix of vertices in linear network representation
  A <- L$m*1
  # vertex indices with degree 2
  deg_v <- as.numeric(table(c(L$from, L$to)))
  ind <- which(deg_v == 2)

  # initialize
  ind2 <- NULL
  i <- 0
  P <- list()

  # remove vertices with degree 2 until none is left
  while(length(ind) > 0){
    i <- i+1

    # find the two adjacent vertices to vertex with degree 2
    adj <- which(A[ind[1], ] == 1)
    P[[i]] <- dplyr::tibble(from = c(adj[1], ind[1]),
                            to = c(ind[1], adj[2]),
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
      P[[i]]$m[k] <- which(L$from == P[[i]]$from[k] & L$to == P[[i]]$to[k] | L$from == P[[i]]$to[k] & L$to == P[[i]]$from[k])
      P[[i]]$length[k] <- L$d[P[[i]]$m[k]]
    }

    A[P[[i]]$from[1], P[[i]]$to[nrow(P[[i]])]] <-
      A[P[[i]]$to[nrow(P[[i]])], P[[i]]$from[1]] <- 1
    # add line segments to the delete vector
    ind2 <- unique(c(ind2, P[[i]]$m))

    # remove vertices from the current vector of vertices with degree 2
    ind <- setdiff(ind, P[[i]]$to[1:(nrow(P[[i]])-1)])
    G$lins$e[P[[i]]$m] <- i
  }
  ind <- which(deg_v == 2)
  G$W <- L$vertices$n - length(ind)
  G$M <-  L$lines$n - length(ind2) + length(P)

  G$vertices$v[setdiff(1:nrow(G$vertices), ind)] <- 1:G$W
  G$lins$e[which(is.na(G$lins$e))] <- (length(P) + 1):G$M

  G$d <- G$lins %>%
    dplyr::group_by(e) %>%
    dplyr::summarize(length = sum(length), .groups = "drop") %>%
    dplyr::pull(length)

  # delete rows and cols in A
  A <- A[-ind, ]
  A <- A[, -ind]
  G$adjacency <- upper.tri(A)*A
  net <- network::as.network(G$adjacency)
  G$incidence <- network::as.matrix.network(net, matrix.type = "incidence")

  class(G) <- "gn"
  G
}

#' Methods for Geometric Networks
#'
#' The function as.gnds converts an object of class gnd to an object of
#' class gnds
#'
#' @param x an object of class gnd or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @return an object of class gnds
#' @export
as.gnds <- function(x, delta = NULL, h = NULL, r = 1){
  if (!inherits(x, c("gn", "gnd", "linnet", "lpp"))){
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
  class(x) <- c(class(x), "gns")
  x
}
