#' Transmute to Geometric Network
#'
#' \code{as_gn} transmutes an existing object into a geometric network, an object
#' of class \code{gn}.
#'
#' @param x An object that could reasonably be transmuted to an
#' object of class \code{gn}.
#' @param ... Other arguments.
#' @return An object of class \code{gn}.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export
#' @examples
#' library(spatstat.data)
#' G <- as_gn(simplenet)
#' summary(G)

as_gn <- function(x, ...){
  UseMethod("as_gn")
}

#' Transmute to Point Pattern on a Geometric Network
#'
#' \code{as_gnpp} transmutes an existing object into a point pattern on a geometric
#' network, an object of class \code{gnpp}.
#'
#' @param x An object that could reasonably transmuted to an object of class
#' \code{gnpp}.
#' @param ... Other arguments.
#' @return An object of class \code{gnpp}.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export
#' @examples
#' library(spatstat.data)
#' X <- as_gnpp(chicago)
#' summary(X)

as_gnpp <- function(x, ...){
  UseMethod("as_gnpp")
}

#' Transmute to Point Pattern on a Linear Network
#'
#' \code{as_lpp} transmutes an existing object into a point pattern on a
#' geometric network, an object of class \code{lpp}.
#'
#' @param x An object that could reasonably transmuted to an object of class
#' \code{lpp}.
#' @param ... Other arguments.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @return A point pattern on a linear network, an object of class \code{lpp}.
#' @export

as_lpp <- function(x, ...){
  UseMethod("as_lpp")
}


#' @param spatstat Set to \code{TRUE} if retransformation to an object of any
#' \code{spatstat} class is desired. Some elements of these objects (such as
#' the window) are otherwise discarded when being transmuted to an object of
#' class \code{gn}.
#' @rdname as_gn
#' @export

as_gn.linnet <- function(x, ..., spatstat = FALSE){
  stopifnot(inherits(spatstat, "logical"))
  l <- a1 <- a2 <- e <- NULL
  L <- x
  d <- diag(L$dpath[L$from, L$to])
  G <- list(
    vertices = tibble(a = 1:L$vertices$n,
                             v = NA,
                             x = L$vertices$x,
                             y = L$vertices$y),
    lins = NULL,
    adjacency = NULL, incidence = NULL,
    d = NULL, q = 2, W = NULL, M = NULL,
    unit = L$window$units
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
    curves[[i]] <- tibble(a1 = c(adj[1], v_deg2[1]),
                     a2 = c(v_deg2[1], adj[2]),
                     l = NA,
                     length = NA)
    # go into the direction of the first adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[1], ]) == 2){
      adj.new <- which(A[adj[1], ] == 1)
      left <- adj.new[which(!is.element(adj.new, curves[[i]]$a1))]
      curves[[i]] <- rbind(c(left, adj[1], NA, NA), curves[[i]])
      adj[1] <- left
    }
    # go into the direction of the second adjecent vertex and search for more
    # vertices with degree 2
    while(sum(A[adj[2], ]) == 2){
      adj.new <- which(A[adj[2], ] == 1)
      right <- adj.new[which(!is.element(adj.new, curves[[i]]$a2))]
      curves[[i]] <- rbind(curves[[i]], c(adj[2], right, NA, NA))
      adj[2] <- right
    }
    # save the line indices and their corresponding lengths which are removed from the network
    for (k in 1:nrow(curves[[i]])) {
      curves[[i]]$l[k] <- which(L$from == curves[[i]]$a1[k] & L$to == curves[[i]]$a2[k] |
                                    L$from == curves[[i]]$a2[k] & L$to == curves[[i]]$a1[k])
      curves[[i]]$length[k] <- d[curves[[i]]$l[k]]
    }
    # add new connection to adjacancy matrix
    A[curves[[i]]$a1[1], curves[[i]]$a2[nrow(curves[[i]])]] <-
      A[curves[[i]]$a2[nrow(curves[[i]])], curves[[i]]$a1[1]] <- 1
    # add line segments to the delete vector
    lins_remove <- unique(c(lins_remove, curves[[i]]$l))
    # remove vertices from the current vector of vertices with degree 2
    v_deg2 <- setdiff(v_deg2, curves[[i]]$a2[1:(nrow(curves[[i]])-1)])
    # add information for geometric network representation
    cs <- cumsum(curves[[i]]$length)
    curves[[i]]$e <- i
    curves[[i]]$frac1 <- c(0, cs[-length(cs)]/cs[length(cs)])
    curves[[i]]$frac2 <- curves[[i]]$length/sum(curves[[i]]$length)
  }
  v_deg2 <- which(degrees == 2)
  G$W <- L$vertices$n - length(v_deg2)
  G$M <-  L$lines$n - length(lins_remove) + length(curves)
  curves <- bind_rows(curves) %>% mutate(a1_x = L$vertices$x[a1],
                                   a1_y = L$vertices$y[a1],
                                   a2_x = L$vertices$x[a2],
                                   a2_y = L$vertices$y[a2])
  ind_lins <- setdiff(1:L$lines$n, lins_remove)
  lins <- tibble(l = ind_lins, e = (i+1):G$M,
                  a1 = L$from[ind_lins], a2 = L$to[ind_lins],
                  a1_x = L$vertices$x[a1], a1_y = L$vertices$y[a1],
                  a2_x = L$vertices$x[a2], a2_y = L$vertices$y[a2],
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
  if (spatstat) {
    G$window <- x$window
    lins <- G$lins %>% arrange(l)
    G$seg_permut <- which(lins$a1_x == L$lines$end$x1)
  }
  class(G) <- "gn"
  G
}

#' @rdname as_gn
#' @export

as_gn.gnpp <- function(x, ...){
  G <- x$network
  class(G) <- "gn"
  G
}

#' @rdname as_gn
#' @export

as_gn.gnppfit <- function(x, ...){
  G <- x$network
  class(G) <- "gn"
  G
}

#' @rdname as_gn
#' @import spatstat.linnet
#' @export

as_gn.lpp <- function(x, ..., spatstat = FALSE) {
  if (!is.logical(spatstat)) stop ("spatstat must be logical")
  L <- as.linnet(x)
  as_gn(L, spatstat = spatstat)
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
#' the window) are otherwise discarded when being transmuted to an object of
#' class \code{gn}.
#' @import dplyr
#' @importFrom spatstat.linnet as.linnet
#' @export

as_gnpp.lpp <- function(x, ..., spatstat = FALSE){
  stopifnot(inherits(spatstat, "logical"))
  frac1 <- tp_e <- tp_l <- frac2 <- e <- y <- l <- xx <-  NULL
  G <- as_gn(as.linnet(x), spatstat = spatstat)
  data <- tibble(l = x$data$seg, tp_l = x$data$tp,
                xx = x$data$x, y = x$data$y)
  if (ncol(x$data) > 4){
    covariates <- as_tibble(as.data.frame(x$data)[, -(1:4)])
    colnames(covariates) <- colnames(x$data)[-(1:4)]
    data <- bind_cols(data, covariates)
    data <- left_join(data, G$lins, by = "l") %>%
      mutate(tp_e = frac1 + tp_l*frac2, x = xx) %>%
      select(l, tp_l, e, tp_e, x, y, colnames(covariates))
  } else {
    data <- left_join(data, G$lins, by = "l") %>%
      mutate(tp_e = frac1 + tp_l*frac2, x = xx) %>%
      select(l, tp_l, e, tp_e, x, y)
  }
  if (!spatstat) data <- data %>% arrange(e, l)
  X <- list(data = data, network = G)
  class(X) <- "gnpp"
  X
}

#' Transmute to Linear Network
#'
#' \code{as.linnet.gn} is a method for the generic function
#' \code{\link[spatstat.linnet]{as.linnet}} which transmutes a geometric network
#' (object of class \code{gn}) to a linear network (object of
#' class \code{\link[spatstat.linnet]{linnet}}).
#'
#' @param X An object of class \code{gn}.
#' @param ... Other arguments.
#' @import spatstat.linnet spatstat.geom
#' @import dplyr
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @return A linear network, an object of class \code{linnet}.
#' @export
#' @examples
#' library(spatstat.data)
#' library(spatstat.linnet)
#'
#' x <- as.linnet(small_gn)
#' plot(x)
#'
#' L <- simplenet
#' X <- as_gn(L, spatstat = TRUE)
#' x <- as.linnet(X)
#' # TRUE
#' all.equal(x, L)

as.linnet.gn <- function(X, ...) {
  l <- a2 <- NULL
  window <- owin(xrange = range(X$vertices$x),
                 yrange = range(X$vertices$y))
  vertices <- ppp(x = X$vertices$x, y = X$vertices$y, window = window)
  edges <- X$lins %>% arrange(l) %>% select(a1, a2) %>% as.matrix()
  if (!is.null(X$window)) {
    vertices$window <- X$window
    a1 <- edges[, 1]
    edges[X$seg_permut, 1] <- edges[X$seg_permut, 2]
    edges[X$seg_permut, 2] <- a1[X$seg_permut]
  }
  out <- linnet(vertices = vertices, edges = edges)
}

#' @rdname as_lpp
#' @param x An object of class \code{gnpp}.
#' @import spatstat.linnet spatstat.geom
#' @import dplyr
#' @export
#' @examples
#' library(spatstat.data)
#' library(spatstat.linnet)
#'
#' x <- as_lpp(montgomery)
#' plot(x)
#'
#' L <- chicago
#' X <- as_gnpp(chicago, spatstat = TRUE)
#' x <- as_lpp(X)
#' # TRUE
#' all.equal(x, L)

as_lpp.gnpp <- function(x, ...) {
  L <- as.linnet(as_gn(x))
  data <- x$data
  marks <- NULL
  if (!is.null(data$marks)) marks <- data$marks
  out <- as.lpp(x = data$x, y = data$y,
                seg = data$l, tp = data$tp_l, L = L, marks = marks)
  out
}

#' @rdname as_gn

as_gn.lppfit <- function(x, ...) {
  x$network
}
