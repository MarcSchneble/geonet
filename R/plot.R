#' Plot Method for Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gn}.
#'
#' Details
#'
#' @param x A geometric network (object of class \code{gn}).
#' @param ... Further arguments passed to plot.
#' @param title Main title of the plot.
#' @param title_x x-axis title (ignored if \code{frame = FALSE}).
#' @param title_y y-axis title (ignored if \code{frame = FALSE}).
#' @param size Plotting size of the line segments. This specificies the
#' \code{size} argument of \code{\link[ggplot2]{geom_segment}}.
#' @param frame Should a frame be drawn around the network?
#' @return Invisibly returns an object of class \code{ggplot}.
#' @import ggplot2
#' @export

plot.gn <- function(x, ..., title = "", title_x = "", title_y = "",
                    size = 1, frame = TRUE) {
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  g <- ggplot(x$lins)
  if (!frame) {
    g <- g + theme_void()
  } else {
    g <- g + theme_bw()
  }
  g <- g + geom_segment(aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = size,
                        lineend = "round", linejoin = "bevel") +
    labs(x = title_x, y = title_y, title = title) +
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
#' @param x A point pattern on a geometric network
#' (object of class \code{gnpp}).
#' @param ... Further arguments passed to plot.
#' @param covariate Character vector of length 1, name of the covariate
#' that should be plotted.
#' @param title Main title of the plot.
#' @param title_x x-axis title (ignored if \code{frame = FALSE}).
#' @param title_y y-axis title (ignored if \code{frame = FALSE}).
#' @param size_lines Plotting size of the line segments. This specificies the
#' \code{size} argument of \code{\link[ggplot2]{geom_segment}}.
#' @param size_points Plotting size of the point pattern. This specificies the
#' \code{size} argument of \code{\link[ggplot2]{geom_point}}.
#' @param frame Should a frame be drawn around the network?
#' @return Invisibly returns an object of class \code{ggplot}.
#' @import ggplot2
#' @export

plot.gnpp <- function(x, ..., covariate = NULL,
                      title = "", title_x = "", title_y = "",
                      size_lines = 1, size_points = 1, frame = TRUE) {
  if (length(covariate) > 1){
    stop("Currently, only one covariate can be plotted.")
  }
  y <- v1_x <- v1_y <- v2_x <- v2_y <- NULL
  g <- ggplot(x$network$lins)
  if (!frame) {
    g <- g + theme_void()
  } else {
    g <- g + theme_bw()
  }
  g <- g + geom_segment(aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y),
                          size = size_lines, lineend = "round", linejoin = "mitre") +
    labs(x = title_x, y = title_y, title = title) +
    theme(panel.grid = element_blank(),
                   plot.title = element_text(hjust = 0.5))
  if (is.null(covariate)) {
    g <- g + geom_point(data = x$data, aes(x = x, y = y),
                        color = "red", size = size_points)
  } else {
    g <- g + geom_point(data = x$data, aes(x = x, y = y, color = !!sym(covariate)),
                        size = size_points)
  }
  print(g)
  invisible(g)
}

#' Plot Method for Fitted Models on Geometric Networks
#'
#' This is the \code{plot} method for an object of class \code{gnppfit}.
#'
#' Details
#'
#' @param x A model fit on a geometric network
#' (an object of class \code{gnppfit}).
#' @param ... Further arguments passed to plot
#' @param select Allows the plot for a single model term to be selected for
#' printing. e.g. if you just want the plot for the second smooth term set
#' \code{select = 2}.
#' @param title Main title of the plot.
#' @param title_x x-axis title (ignored if \code{frame = FALSE}).
#' @param title_y y-axis title (ignored if \code{frame = FALSE}).
#' @param size Plotting size of the line segments. This specificies the
#' \code{size} argument of \code{\link[ggplot2]{geom_segment}}.
#' @param frame Should a frame be drawn around the network?
#' @param sol asd
#' @return Invisibly returns a list of objects with elements of class ggplot.
#' @import ggplot2
#' @importFrom  grDevices devAskNewPage
#' @export

plot.gnppfit <- function(x, ..., select = NULL, title = "", title_x = "x", title_y = "y",
                         size = 1, frame = TRUE, sol = 100) {

  e <- xend <- y <- yend <- intensity <- lwr <- upr <- NULL

  stopifnot(inherits(x, "gnppfit"))

  df <- tibble(seg = integer(0), e = integer(0),
                      x = numeric(0), xend = numeric(0),
                      y = numeric(0), yend = numeric(0),
                      z = numeric(0))
  G <- as.gn(x)

  g <- df <- vector("list", 1 + length(x$effects$smooth))

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
      df[[1]] <- bind_rows(df, tibble(seg = dat$seg[i], e = m,
                                               x = utils::head(xx, -1), xend = xx[-1],
                                               y = utils::head(yy, -1), yend = yy[-1],
                                               z = zz))
    }
  }

  # get intensity
  B <- getBplot(x, df[[1]])
  ind <- grep("G.", x$names_theta)
  df[[1]]$intensity <- as.vector(exp(B[, ind]%*%x$theta[ind]))
  #

  g[[1]] <- ggplot(df[[1]])
  if (!frame) {
    g[[1]] <- g[[1]] + theme_void()
  } else {
    g[[1]] <- g[[1]] + theme_bw()
  }
  g[[1]] <- g[[1]] + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = intensity),
                          size = size, lineend = "round", linejoin = "mitre") +
    labs(color = "Intensity", x = title_x, y = title_y, title = title) +
    theme(panel.grid = element_blank(),
                 plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red")

  if (length(g) > 1){
    for (i in 2:length(g)) {
      df[[i]] <- x$effects$smooth[[i - 1]]
      g[[i]] <- ggplot(df[[i]]) +
        geom_ribbon(aes(x = x, ymin = exp(lwr), ymax = exp(upr)), fill = "grey50") +
        geom_line(aes(x = x, y = exp(y)), color = "red") +
        theme_bw()
    }
  }
  ask <- TRUE
  for (i in 1:length(g)) {
    print(g[[i]])
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      ask <- FALSE
    }
  }
  invisible(g)
}
