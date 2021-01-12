#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
as.geonet <- function(object){
  if (!inherits(object, "linnet")){
    object <- spatstat::as.linnet(object)
    if (!inherits(object, "linnet")){
      stop("Object must be of class 'linnet' or must be converable to an object of class 'linnet'")
    }
  }
  net <- list()
  net$vertices <- tibble(x = object$vertices$x, y = object$vertices$y)

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
    P[[i]] <- tibble(from = c(adj[1], ind.deg.2[1]),
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
    curves[[k]] <- tibble(from = object$from[m], to = object$to[m], m = k, length = object$d[m])
  }
  net$curves <- bind_rows(curves)

  class(net) <- "geonet"
  net
}
