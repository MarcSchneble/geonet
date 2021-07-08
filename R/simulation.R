#' Simulate Uniform Points on a Geometric Network
#'
#' \code{runifgn} simulates points on a geometric network according
#' to a uniform density.
#'
#' @param n Number of random points to generate. A nonnegative integer.
#' @param G A geometric network (object of class \code{gn}).
#' @return An object of class \code{gnpp}.
#' @import dplyr
#' @importFrom stats runif
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export


runifgn <- function(n, G){
  e <- seg <- tp <- frac1 <- frac2 <- x <- y <- NULL
  lins <- G$lins %>% arrange(id)
  cumlen <- c(0, cumsum(lins$length))
  dx <- lins$v2_x - lins$v1_x
  dy <- lins$v2_y - lins$v1_y
  uu <- runif(n, min = 0, max = sum(G$d))
  ii <- findInterval(uu, cumlen, rightmost.closed = TRUE, all.inside = TRUE)
  tt <- (uu - cumlen[ii])/lins$length[ii]
  xx <- lins$v1_x[ii] + tt*dx[ii]
  yy <- lins$v1_y[ii] + tt*dy[ii]
  data <- tibble(id = ii, tp = tt)
  data <- left_join(data, lins, by = "id") %>%
    mutate(x = xx, y = yy, tp = frac1 + tp*frac2) %>%
    select(id, e, tp, x, y) %>%
    arrange(e, tp)
  out <- list(data = data, network = G)
  class(out) <- "gnpp"
  out
}
