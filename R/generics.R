#' Methods for Geometric Networks
#'
#' This is the print method for object of the class gn.
#'
#' @param x an object of class gn
#' @param ... further arguments passed to print
#' @export
print.gn <- function(x, ...){
  stopifnot(inherits(x, "gn"))
  cat(paste("Geometric network in", x$q, "dimensions"))
  invisible(NULL)
}

#' Methods for Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gn}.
#'
#' Details
#'
#' @param x an object of class \code{gn}
#' @param ... further arguments passed to plot
#' @param x.title x-axis title
#' @param y.title y-axis title
#' @param title plot title
#' @return an object of class ggplot
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @export
plot.gn <- function(x, ..., x.title = "x", y.title = "y", title = NULL){
  from <- to <- NULL
  x.from <- y.from <- x.to <- y.to <- NULL
  dat.lins <- x$lins %>%
    dplyr::mutate(x.from = x$vertices$x[from],
                  y.from = x$vertices$y[from],
                  x.to = x$vertices$x[to],
                  y.to = x$vertices$y[to])
  g <- ggplot2::ggplot(dat.lins) +
    ggplot2::geom_segment(aes(x = x.from, y = y.from, xend = x.to, yend = y.to)) +
    ggplot2::labs(x = x.title, y = y.title, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = element_blank(),
                   plot.title = element_text(hjust = 0.5))
  print(g)
  invisible(g)
}
