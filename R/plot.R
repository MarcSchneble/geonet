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
#' @param size asd
#' @return an object of class ggplot
#' @import ggplot2
#' @export
plot.gn <- function(x, ..., x.title = "x", y.title = "y", title = NULL, size = 1){
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  g <- ggplot(x$lins) +
    geom_segment(aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = size,
                          lineend = "round", linejoin = "bevel") +
    labs(x = x.title, y = y.title, title = title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
                   plot.title = element_text(hjust = 0.5))
  print(g)
  invisible(g)
}

#' Plot Method for Points Patterns on Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gnpp}.
#'
#' Details
#'
#' @param x A geometric network with data (an object of class \code{gnpp}).
#' @param ... further arguments passed to plot
#' @param title plot title
#' @param title_y asd
#' @param title_x asd
#' @param size asdas
#' @return an object of class ggplot
#' @import ggplot2
#' @export

plot.gnpp <- function(x, ..., title_x = "x", title_y = "y", title = "", size = 1){
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  y <- NULL
  g <- ggplot(x$lins) +
    geom_segment(aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y),
                          size = size, lineend = "round", linejoin = "mitre") +
    labs(x = title_x, y = title_y, title = title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
                   plot.title = element_text(hjust = 0.5)) +
    geom_point(data = x$data, aes(x = x, y = y), color = "red")
  print(g)
  invisible(g)
}

#' Plot Method for Points Patterns on Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gnpp}.
#'
#' Details
#'
#' @param x A geometric network with data (an object of class \code{gnppfit}).
#' @param ... further arguments passed to plot
#' @param title plot title
#' @param title_y asd
#' @param title_x asd
#' @param size asdas
#' @param sol asd
#' @return an object of class ggplot
#' @import ggplot2
#' @export

plot.gnppfit <- function(x, ..., title_x = "x", title_y = "y", title = "", size = 1, sol = 100) {

  e <- xend <- y <- yend <- intensity <- NULL

  stopifnot(inherits(x, "gnppfit"))

  df <- tibble(seg = integer(0), e = integer(0),
                      x = numeric(0), xend = numeric(0),
                      y = numeric(0), yend = numeric(0),
                      z = numeric(0))
  G <- as.gn(x)

  for (m in 1:G$M) {
    dat <- filter(G$lins, e == m)
    cs <- c(0, cumsum(dat$length))
    dx <- dat$v2_x - dat$v1_x
    dy <- dat$v2_y - dat$v1_y
    for (i in 1:length(dat$seg)) {
      tt <- seq(0, 1, 1/sol)
      xx <- dat$v1_x[i] + tt*dx[i]
      yy <- dat$v1_y[i] + tt*dy[i]
      zz <- cs[i] + (tt - 1/(2*sol))[-1]*dat$length[i]
      df <- bind_rows(df, tibble(seg = dat$seg[i], e = m,
                                               x = utils::head(xx, -1), xend = xx[-1],
                                               y = utils::head(yy, -1), yend = yy[-1],
                                               z = zz))
    }
  }

  # get intensity
  B <- getBplot(x, df)
  ind <- grep("G.", x$names_theta)
  df$intensity <- as.vector(exp(B[, ind]%*%x$theta[ind]))
  #

  g <- ggplot(df) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = intensity),
                          size = size, lineend = "round", linejoin = "mitre") +
    labs(color = "Intensity", x = title_x, y = title_y, title = title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
                 plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red")
  print(g)
  invisible(df)
}
