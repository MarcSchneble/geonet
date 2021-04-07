#' Coerce to Geometric Network
#'
#' \code{as_gn} coerces an existing object into a geometric network, an object
#' of class \code{gn}.
#'
#' @param x An object that could reasonably be coerced to a geometric network
#' (object of class gn).
#' @param ... Further arguments passed to \code{as_gn}.
#' @export

as_gn <- function(x, ...){
  UseMethod("as_gn")
}

#' Coerece to Point Pattern on a Geometric Network
#'
#' \code{as_gnpp} coerces an existing object into a point pattern on a geometric
#' network, an object
#' of class \code{gnpp}.
#'
#' @param x An object that could reasonable coerced to a point pattern on a
#' geometric network (object of class \code{gnpp}).
#' @param ... further
#' @export

as_gnpp <- function(x, ...){
  UseMethod("as_gnpp")
}



#' @param spatstat Set to \code{TRUE} if retransformation to an object of any
#' \code{spatstat} class is desired. Some elements of these objects (such as
#' the window) are discared when being coerced to an object of class
#' \code{gn}.
#' @rdname as_gn
#' @export

as_gn.linnet <- function(x, ..., spatstat = FALSE){
  if (!inherits(x, "linnet")) stop("x muss be of class 'linnet'")
  v1 <- v2 <- e <- NULL
  L <- x
  d <- diag(L$dpath[L$from, L$to])
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
  degrees <- as.numeric(table(c(L$from, L$to)))
  v_deg2 <- which(degrees == 2)
  # initialize
  lins_remove <- NULL
  curves <- list()
  i <- 0
  # remove vertices with degree 2 until none is left
  while(length(v_deg2) > 0){
    i <- i+1
    # find the two adjacent vertices to vertex with degree 2
    adj <- which(A[v_deg2[1], ] == 1)
    curves[[i]] <- tibble(v1 = c(adj[1], v_deg2[1]),
                     v2 = c(v_deg2[1], adj[2]),
                     seg = NA,
                     length = NA)
    # go into the direction of the first adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[1], ]) == 2){
      adj.new <- which(A[adj[1], ] == 1)
      l <- adj.new[which(!is.element(adj.new, curves[[i]]$v1))]
      curves[[i]] <- rbind(c(l, adj[1], NA, NA), curves[[i]])
      adj[1] <- l
    }
    # go into the direction of the second adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[2], ]) == 2){
      adj.new <- which(A[adj[2], ] == 1)
      r <- adj.new[which(!is.element(adj.new, curves[[i]]$v2))]
      curves[[i]] <- rbind(curves[[i]], c(adj[2], r, NA, NA))
      adj[2] <- r
    }
    # save the line indices and their corresponding lengths which are removed from the network
    for (k in 1:nrow(curves[[i]])) {
      curves[[i]]$seg[k] <- which(L$from == curves[[i]]$v1[k] & L$to == curves[[i]]$v2[k] |
                                    L$from == curves[[i]]$v2[k] & L$to == curves[[i]]$v1[k])
      curves[[i]]$length[k] <- d[curves[[i]]$seg[k]]
    }
    # add new connection to adjacancy matrix
    A[curves[[i]]$v1[1], curves[[i]]$v2[nrow(curves[[i]])]] <-
      A[curves[[i]]$v2[nrow(curves[[i]])], curves[[i]]$v1[1]] <- 1
    # add line segments to the delete vector
    lins_remove <- unique(c(lins_remove, curves[[i]]$seg))
    # remove vertices from the current vector of vertices with degree 2
    v_deg2 <- setdiff(v_deg2, curves[[i]]$v2[1:(nrow(curves[[i]])-1)])
    # add information for geometric network representation
    cs <- cumsum(curves[[i]]$length)
    curves[[i]]$e <- i
    curves[[i]]$frac1 <- c(0, cs[-length(cs)]/cs[length(cs)])
    curves[[i]]$frac2 <- curves[[i]]$length/sum(curves[[i]]$length)
  }
  v_deg2 <- which(degrees == 2)
  G$W <- L$vertices$n - length(v_deg2)
  G$M <-  L$lines$n - length(lins_remove) + length(curves)
  curves <- bind_rows(curves) %>% mutate(v1_x = L$vertices$x[v1],
                                   v1_y = L$vertices$y[v1],
                                   v2_x = L$vertices$x[v2],
                                   v2_y = L$vertices$y[v2])
  ind_lins <- setdiff(1:L$lines$n, lins_remove)
  lins <- tibble(seg = ind_lins, e = (i+1):G$M,
                  v1 = L$from[ind_lins], v2 = L$to[ind_lins],
                  v1_x = L$vertices$x[v1], v1_y = L$vertices$y[v1],
                  v2_x = L$vertices$x[v2], v2_y = L$vertices$y[v2],
                  length = diag(L$dpath[L$from, L$to])[ind_lins],
                  frac1 = 0, frac2 = 1)
  G$lins <- bind_rows(lins, curves)
  G$lins <- G$lins %>% arrange(e)
  G$vertices$v[setdiff(1:nrow(G$vertices), v_deg2)] <- 1:G$W
  G$d <- G$lins %>% group_by(e) %>%
    summarize(length = sum(length), .groups = "drop") %>% pull(length)
  # delete rows and cols in A which relate to vertices with degree 2
  A <- A[-v_deg2, ]
  G$adjacency <- A[, -v_deg2]
  G$incidence <- incidence(G$vertices, G$lins)
  class(G) <- "gn"
  G
}

#' @rdname as_gn
#' @export

as_gn.gnpp <- function(x, ...){
  if (!inherits(x, "gnpp")){
    stop("Object must be of class 'gpp'")
  }
  G <- x$network
  class(G) <- "gn"
  G
}

#' @rdname as_gn
#' @export

as_gn.gnppfit <- function(x, ...){
  if (!inherits(x, "gnppfit")){
    stop("Object must be of class 'gnppfit'")
  }
  G <- x$network
  class(G) <- "gn"
  G
}

#' @rdname as_gnpp
#' @export

as_gnpp.gnppfit <- function(x, ...){
  X <- list(data = x$data, network = x$network)
  class(X) <- "gnpp"
  X
}


#' @rdname as_gnpp
#' @param spatstat Set to \code{TRUE} if retransformation to an object of any
#' \code{spatstat} class is desired. Some elements of these objects (such as
#' the window) are discared when being coerced to an object of class
#' \code{gn}.
#' @import dplyr
#' @importFrom spatstat.linnet as.linnet
#' @export

as_gnpp.lpp <- function(x, ..., spatstat = FALSE){
  frac1 <- tp <- frac2 <- e <- y <- seg <- xx <-  NULL
  if (!inherits(x, "lpp")){
      stop("Object must be of class 'lpp'")
  }
  G <- as_gn(as.linnet(x))
  data <- tibble(seg = x$data$seg, tp = x$data$tp,
                xx = x$data$x, y = x$data$y)
  if (ncol(x$data) > 4){
    covariates <- as_tibble(as.data.frame(x$data)[, -(1:4)])
    colnames(covariates) <- colnames(x$data)[-(1:4)]
    data <- bind_cols(data, covariates)
  }
  data <- left_join(data, G$lins, by = "seg") %>%
    mutate(tp = frac1 + tp*frac2, x = xx) %>%
    select(seg, e, tp, x, y, colnames(covariates)) %>%
    arrange(e, tp)
  X <- list(data = data, network = G)
  class(X) <- "gnpp"
  X
}

as.linnet.gn <- function(X, ...) {

}
