#' Simluate Uniform Points on a Geometric Network
#'
#' The function runifgn simulates points on a geometric network
#'
#' @param n Number of random points to generate. A nonnegative integer.
#' @param G A geometric network (object of class "gn", see gn)
#' @return an object of class gnd
#' @export
#' @importFrom dplyr %>%

runifgn <- function(n, G){
  e <- id <- tp <- NULL
  cumlen <- c(0, cumsum(G$lins$length))
  dx <- G$lins$v2_x - G$lins$v1_x
  dy <- G$lins$v2_y - G$lins$v1_y
  uu <- stats::runif(n, min = 0, max = sum(G$d))
  ii <- findInterval(uu, cumlen, rightmost.closed = TRUE, all.inside = TRUE)
  tt <- (uu - cumlen[ii])/G$lins$length[ii]
  xx <- G$lins$v1_x[ii] + tt*dx[ii]
  yy <- G$lins$v1_y[ii] + tt*dy[ii]
  dat <- dplyr::tibble(id = ii, tp = tt)
  dat <- dplyr::left_join(dat, G$lins, by = "id") %>%
    dplyr::mutate(x = xx, y = yy) %>%
    dplyr::select(id, e, tp, x, y) %>%
    dplyr::arrange(e, tp, id)
  G$data <- dat
  class(G) <- "gnd"
  G
}
