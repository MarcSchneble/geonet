#' Generic function
#'
#' Gasdda
#'
#' @param x An object that should be converted to a geometric network
#' (object of class gn)

as.gn <- function(x, ...){
  UseMethod("as.gn")
}

#' Generic function
#'
#' Gasdda
#'
#' @param x An object that should be converted to a point pattern on a
#' geometric network (object of class gnpp).

as.gnpp <- function(x, ...){
  UseMethod("as.gnpp")
}

#' Methods for Linear Networks
#'
#' The function as.gn converts an object of class linnet to an object of
#' class gn
#'
#' @param object an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @return an object of class gn
#' @export
as.gn.linnet <- function(x){
  if (!inherits(x, "linnet")){
    x <- spatstat::as.linnet(x)
    if (!inherits(x, "linnet")){
      stop("Object must be of class 'linnet' or must be converable to an object of class 'linnet'")
    }
  }
  L <- x
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
                         length = diag(L$dpath[L$from, L$to]),
                         frac1 = 0,
                         frac2 = 1),
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

    # get fraction of line segments with respect to the whole curce at the end and at the beginning
    cs <- cumsum(P[[i]]$length)
    G$lins$frac1[which(G$lins$e == i)] <- c(0, cs[-length(cs)]/cs[length(cs)])
    G$lins$frac2[which(G$lins$e == i)] <- cs/cs[length(cs)]
  }
  ind <- which(deg_v == 2)
  G$W <- L$vertices$n - length(ind)
  G$M <-  L$lines$n - length(ind2) + length(P)

  G$vertices$v[setdiff(1:nrow(G$vertices), ind)] <- 1:G$W
  G$lins$e[which(is.na(G$lins$e))] <- (length(P) + 1):G$M
  G$lins <- G$lins %>% dplyr::arrange(e, id)

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

#' Methods for Point Patterns on Linear Networks
#'
#' The function as.gnpp.lpp converts an object of class linnet to an object of
#' class gn
#'
#' @param object an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @return an object of class gn
#' @export

as.gnpp.lpp <- function(x, ...){
  if (!inherits(x, "lpp")){
    x <- spatstat::as.lpp(x)
    if (!inherits(x, "lpp")){
      stop("Object must be of class 'lpp' or must be converable to an object of class 'lpp'")
    }
  }
  G <- as.gn(spatstat::as.linnet(x))
  X <- x
  dat <- dplyr::tibble(id = X$data$seg,
                       tp = X$data$tp,
                       x = X$data$x,
                       y = X$data$y)
  G$data <- dplyr::left_join(dat, G$lins, by = "id") %>%
    dplyr::select(id, e, tp, x, y) %>%
    dplyr::arrange(e, id, tp)
  class(G) <- "gnpp"
  G
}

