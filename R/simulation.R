#' Simluate Uniform Points on a Geometric Network
#'
#' The function runifgn simulates points on a geometric network
#'
#' @param n Number of random points to generate. A nonnegative integer.
#' @param G A geometric network (object of class "gn", see gn)
#' @return an object of class gnd
#' @import dplyr
#' @importFrom stats runif
#' @export


runifgn <- function(n, G){
  e <- seg <- tp <- frac1 <- frac2 <- x <- y <- NULL
  cumlen <- c(0, cumsum(G$lins$length))
  dx <- G$lins$v2_x - G$lins$v1_x
  dy <- G$lins$v2_y - G$lins$v1_y
  uu <- runif(n, min = 0, max = sum(G$d))
  ii <- findInterval(uu, cumlen, rightmost.closed = TRUE, all.inside = TRUE)
  tt <- (uu - cumlen[ii])/G$lins$length[ii]
  xx <- G$lins$v1_x[ii] + tt*dx[ii]
  yy <- G$lins$v1_y[ii] + tt*dy[ii]
  data <- tibble(id = ii, tp = tt)
  data <- left_join(data, G$lins, by = "id") %>%
    mutate(x = xx, y = yy, tp = frac1 + tp*frac2) %>%
    select(id, e, tp, x, y) %>%
    arrange(e, tp)
  out <- list(data = data, network = G)
  class(out) <- "gnpp"
  out
}
