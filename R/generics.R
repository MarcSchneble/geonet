#' Print Method for Geometric Networks
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

#' Print Method for Geometric Networks Splines
#'
#' This is the print method for object of the class gns.
#'
#' @param x an object of class gns
#' @param ... further arguments passed to print
#' @export
print.gns <- function(x, ...){
  stopifnot(inherits(x, "gns"))
  cat(paste("Geometric network in", x$q, "dimensions"))
  invisible(NULL)
}

#' Plot Method for Geometric Networks
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
#' @export
plot.gn <- function(x, ..., x.title = "x", y.title = "y", title = NULL){
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  g <- ggplot2::ggplot(x$lins) +
    ggplot2::geom_segment(ggplot2::aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y)) +
    ggplot2::labs(x = x.title, y = y.title, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  print(g)
  invisible(g)
}

#' Plot Method for Data on Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gns}.
#'
#' Details
#'
#' @param x A geometric network with data (an object of class \code{gns}).
#' @param ... further arguments passed to plot
#' @param title plot title
#' @return an object of class ggplot
#' @export

plot.gnd <- function(x, ..., title_x = "x", title_y = "y", title = ""){
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  y <- NULL
  g <- ggplot2::ggplot(x$lins) +
    ggplot2::geom_segment(ggplot2::aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y)) +
    ggplot2::labs(x = title_x, y = title_y, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_point(data = x$data, ggplot2::aes(x = x, y = y))
  print(g)
  invisible(g)
}
