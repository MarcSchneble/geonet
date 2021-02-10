#' Coerce to Geometric Network
#'
#' Gasdda
#'
#' @param x An object that should be converted to a geometric network
#' (object of class gn)
#' @param ... further
#' @export

as.gn <- function(x, ...){
  UseMethod("as.gn")
}

#' Coerece to Point Pattern on Geometric Network
#'
#' Gasdda
#'
#' @param x An object that should be converted to a point pattern on a
#' geometric network (object of class gnpp).
#' @param ... further
#' @export

as.gnpp <- function(x, ...){
  UseMethod("as.gnpp")
}

#' Methods for Linear Networks
#'
#' The function as.gn converts an object of class linnet to an object of
#' class gn
#'
#' @param x an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @param ... further
#' @return A geometric network (object of class \code{gnpp})
#' @import dplyr
#' @importFrom spatstat as.linnet
#' @export
#'
as.gn.linnet <- function(x, ...){
  v1 <- v2 <- e <- NULL
  if (!inherits(x, "linnet")){
    x <- as.linnet(x)
    if (!inherits(x, "linnet")){
      stop("Object must be of class 'linnet' or must be converable to an object of class 'linnet'")
    }
  }
  L <- x
  L$d <- diag(L$dpath[L$from, L$to])
  G <- list(
    vertices = tibble(id = 1:L$vertices$n,
                             v = NA,
                             x = L$vertices$x,
                             y = L$vertices$y),
    lins = NULL,
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
    P[[i]] <- tibble(v1 = c(adj[1], ind[1]),
                     v2 = c(ind[1], adj[2]),
                     seg = NA,
                     length = NA)

    # go into the direction of the first adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[1], ]) == 2){
      adj.new <- which(A[adj[1], ] == 1)
      l <- adj.new[which(!is.element(adj.new, P[[i]]$v1))]
      P[[i]] <- rbind(c(l, adj[1], NA, NA), P[[i]])
      adj[1] <- l
    }

    # go into the direction of the second adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[2], ]) == 2){
      adj.new <- which(A[adj[2], ] == 1)
      r <- adj.new[which(!is.element(adj.new, P[[i]]$v2))]
      P[[i]] <- rbind(P[[i]], c(adj[2], r, NA, NA))
      adj[2] <- r
    }

    # save the line indices and their corresponding lengths which are removed from the network
    for (k in 1:nrow(P[[i]])) {
      P[[i]]$seg[k] <- which(L$from == P[[i]]$v1[k] & L$to == P[[i]]$v2[k] | L$from == P[[i]]$v2[k] & L$to == P[[i]]$v1[k])
      P[[i]]$length[k] <- L$d[P[[i]]$seg[k]]
    }

    A[P[[i]]$v1[1], P[[i]]$v2[nrow(P[[i]])]] <-
      A[P[[i]]$v2[nrow(P[[i]])], P[[i]]$v1[1]] <- 1
    # add line segments to the delete vector
    ind2 <- unique(c(ind2, P[[i]]$seg))

    # remove vertices from the current vector of vertices with degree 2
    ind <- setdiff(ind, P[[i]]$v2[1:(nrow(P[[i]])-1)])

    #
    P[[i]]$e <- i
    cs <- cumsum(P[[i]]$length)
    P[[i]]$frac1 <- c(0, cs[-length(cs)]/cs[length(cs)])
    P[[i]]$frac2 <- P[[i]]$length/sum(P[[i]]$length)

    #G$lins$e[P[[i]]$m] <- i

    # get fraction of line segments with respect to the whole curce at the end and at the beginning

    #G$lins$frac1[which(G$lins$e == i)] <- c(0, cs[-length(cs)]/cs[length(cs)])
    #G$lins$frac2[which(G$lins$e == i)] <- P[[i]]$length/sum(P[[i]]$length)
  }
  ind_v <- which(deg_v == 2)
  G$W <- L$vertices$n - length(ind_v)
  G$M <-  L$lines$n - length(ind2) + length(P)

  lins1 <- bind_rows(P) %>% mutate(v1_x = L$vertices$x[v1],
                                   v1_y = L$vertices$y[v1],
                                   v2_x = L$vertices$x[v2],
                                   v2_y = L$vertices$y[v2])
  ind3 <- setdiff(1:L$lines$n, ind2)
  lins2 <- tibble(seg = ind3,
                  e = (length(P) + 1):G$M,
                  v1 = L$from[ind3],
                  v2 = L$to[ind3],
                  v1_x = L$vertices$x[v1],
                  v1_y = L$vertices$y[v1],
                  v2_x = L$vertices$x[v2],
                  v2_y = L$vertices$y[v2],
                  length = diag(L$dpath[L$from, L$to])[ind3],
                  frac1 = 0, frac2 = 1)
  G$lins <- bind_rows(lins2, lins1)
  G$lins <- G$lins %>% dplyr::arrange(e)

  G$vertices$v[setdiff(1:nrow(G$vertices), ind_v)] <- 1:G$W


  G$d <- G$lins %>%
    group_by(e) %>%
    summarize(length = sum(length), .groups = "drop") %>%
    pull(length)

  # delete rows and cols in A
  A <- A[-ind_v, ]
  A <- A[, -ind_v]
  #G$adjacency <- A
  G$incidence <- getIncidence(G)
  G$adjacency <- upper.tri(A)*A
  #net <- network::as.network(G$adjacency)
  #G$incidence <- network::as.matrix.network(net, matrix.type = "incidence")

  class(G) <- "gn"
  G
}

#' Methods for Point Patterns on Linear Networks
#'
#' The function as.gnpp.lpp converts an object of class linnet to an object of
#' class gn
#'
#' @param x an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @param ... asdd
#' @return Point pattern on a geometric network (object of class \code{gnpp}).
#' @export

as.gn.gnpp <- function(x, ...){
  if (!inherits(x, "gnpp")){
    stop("Object must be of class 'gpp'")
  }
  G <- x$network
  class(G) <- "gn"
  G
}

#' Methods for Point Patterns on Linear Networks
#'
#' The function as.gnpp.lpp converts an object of class linnet to an object of
#' class gn
#'
#' @param x an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @param ... asdd
#' @return Point pattern on a geometric network (object of class \code{gnpp}).
#' @export

as.gn.gnppfit <- function(x, ...){
  if (!inherits(x, "gnppfit")){
    stop("Object must be of class 'gnppfit'")
  }
  G <- x$network
  class(G) <- "gn"
  G
}


#' Methods for Point Patterns on Linear Networks
#'
#' The function as.gnpp.lpp converts an object of class linnet to an object of
#' class gn
#'
#' @param x an object of the spatstat class linnet or an object that can
#' be converted to an instance of this class
#' @param ... asdd
#' @return Point pattern on a geometric network (object of class \code{gnpp}).
#' @import dplyr
#' @importFrom spatstat as.linnet
#' @export

as.gnpp.lpp <- function(x, ...){
  frac1 <- tp <- frac2 <- e <- y <- seg <- xx <-  NULL
  if (!inherits(x, "lpp")){
      stop("Object must be of class 'lpp'")
  }
  G <- as.gn(as.linnet(x))
  dat <- tibble(seg = x$data$seg, tp = x$data$tp,
                xx = x$data$x, y = x$data$y)
  if (ncol(x$data) > 4){
    covariates <- as_tibble(as.data.frame(x$data)[, -(1:4)])
    colnames(covariates) <- colnames(x$data)[-(1:4)]
    dat <- bind_cols(dat, covariates)
  }
  dat <- left_join(dat, G$lins, by = "seg") %>%
    mutate(tp = frac1 + tp*frac2, x = xx) %>%
    select(seg, e, tp, x, y, colnames(covariates)) %>%
    arrange(e, tp)
  X <- list(data = dat, network = G)
  class(X) <- "gnpp"
  X
}

