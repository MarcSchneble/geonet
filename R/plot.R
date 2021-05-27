#' Plot Methods for Geometric Network related Objects
#'
#' \code{plot} method for geometric networks, point patterns on geometric
#' networks, or a fitted point process.
#'
#' @param x An object which is related to a geometric network
#' (object of class \code{gn}, \code{gnpp} or \code{gnppfit}).
#' @param ... Other arguments.
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
#' @param sol asd
#' @param scale The scale on which smooth terms should be plotted, either on
#' the log scale (\code{scale = log}, default) or on the exp-scale
#' (\code{scale = exp}).
#' @param select Allows the plot for a single model term to be selected for
#' printing. e.g. if you just want the plot for the second smooth term set
#' \code{select = 2}.
#' @return Invisibly returns an object of class \code{ggplot} or a list
#' of \code{ggplot} objects.
#' @import ggplot2
#' @export

plot.gn <- function(x, ..., title = "", title_x = "", title_y = "",
                    size_lines = 1, frame = TRUE) {
  v1_x <- v1_y <- v2_x <- v2_y <- NULL
  g <- ggplot(x$lins)
  if (!frame) {
    g <- g + theme_void()
  } else {
    g <- g + theme_bw()
  }
  g <- g + geom_segment(aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = size_lines,
                        lineend = "round", linejoin = "bevel") +
    labs(x = title_x, y = title_y, title = title) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  print(g)
  invisible(g)
}

#' @rdname plot.gn
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

#' @rdname plot.gn
#' @import ggplot2
#' @importFrom  grDevices devAskNewPage
#' @importFrom splines splineDesign
#' @export

plot.gnppfit <- function(x, ..., select = NULL, title = "", title_x = "x", title_y = "y",
                         size_lines = 1, frame = TRUE, sol = 100, scale = "log") {

  stopifnot(inherits(x, "gnppfit"))
  e <- xend <- y <- yend <- intensity <- lower <- upper <- NULL
  G <- as_gn(x)
  g <- df <- vector("list", length(x$smooth) + 1)

  for (m in 1:G$M) {
    lins_m <- filter(G$lins, e == m)
    cs <- c(0, cumsum(lins_m$length))
    dx <- lins_m$v2_x - lins_m$v1_x
    dy <- lins_m$v2_y - lins_m$v1_y
    for (i in 1:length(lins_m$id)) {
      tt <- seq(0, 1, 1/sol)
      xx <- lins_m$v1_x[i] + tt*dx[i]
      yy <- lins_m$v1_y[i] + tt*dy[i]
      zz <- cs[i] + (tt - 1/(2*sol))[-1]*lins_m$length[i]
      df[[1]] <- bind_rows(df, tibble(id = lins_m$id[i], e = m,
                                               x = utils::head(xx, -1), xend = xx[-1],
                                               y = utils::head(yy, -1), yend = yy[-1],
                                               z = zz))
    }
  }

  # get baseline intensity of the point pattern on the geometric network
  B <- bspline_design_plot(x, df[[1]])
  theta <- x$coefficients[x$ind[["G"]]]
  df[[1]]$intensity <- as.vector(exp(B%*%theta))
  #

  g[[1]] <- ggplot(df[[1]])
  if (!frame) {
    g[[1]] <- g[[1]] + theme_void()
  } else {
    g[[1]] <- g[[1]] + theme_bw()
  }
  g[[1]] <- g[[1]] + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = intensity),
                          size = size_lines, lineend = "round", linejoin = "mitre") +
    labs(color = "Intensity", x = title_x, y = title_y, title = title) +
    theme(panel.grid = element_blank(),
                 plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red")

  if (length(g) > 1){
    for (i in 2:length(g)) {
      var <- names(x$ind)[i]
      xx <- seq(x$smooth[[var]]$range[1], x$smooth[[var]]$range[2], length.out = 100)
      X <- splineDesign(knots = x$smooth[[var]]$knots, x = xx,
                        ord = x$smooth[[var]]$l + 1, outer.ok = TRUE)
      X <- sweep(X, 2, x$smooth[[var]]$ident)[, -1]
      theta <- x$coefficients[x$ind[[var]]]
      V <- as.matrix(x$V[x$ind[[var]], x$ind[[var]]])
      y <- as.vector(X%*%theta)
      limits <- confidence_band(theta, V, X, R = 1000)
      df[[i]] <- tibble(x = xx, y = y, lower = limits$lower, upper = limits$upper)
      if (scale == "exp"){
        df[[i]] <- mutate(df[[i]], y = exp(y),
                          lower = exp(lower), upper = exp(upper))
      }
      g[[i]] <- ggplot(df[[i]]) +
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "grey50") +
        geom_line(aes(x = x, y = y), color = "red") +
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
