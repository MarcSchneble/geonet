#' Plot Methods for Geometric Network related Objects
#'
#' \code{plot} method for geometric networks, point patterns on geometric
#' networks, or a fitted point process.
#'
#' @param x An object which is related to a geometric network
#' (object of class \code{gn}, \code{gnpp} or \code{gnppfit}).
#' @param ... Other arguments.
#' @param title A named list with names "x", "y" and "plot" which specify the
#' arguments \code{x}, \code{y} and \code{title} of the
#' \code{\link[ggplot2]{labs}} function. Each list entry must
#' be a character vector which has length equal to the number of plots.
#' The list entries can remain unspecified in which case the respective titles
#' are left blank.
#' @param size A named list with names "lines" and "points" which specify the
#' \code{size} argument of \code{\link[ggplot2]{geom_segment}} and
#' \code{\link[ggplot2]{geom_point}}, respectively. Each list entry must be
#' a numeric vector of length one with positive values. If the whole list or one
#' entry remains unspecified, default values are used.
#' @param color A named list with names "lines" and "points" which specify the
#' \code{color} argument of \code{\link[ggplot2]{geom_segment}} and
#' \code{\link[ggplot2]{geom_point}}. Each list entry must specify a valid color.
#' By default, lines and points are plotted in black.
#' @param shape The shape used for plotting the points. An integer between 0 and
#' 25. Default to \code{shape = 1} which shows the points as a circle.
#' @param frame If set to \code{TRUE}, draws a frame around the network and adds
#' tick marks and axis labeling.
#' @param data Set to \code{TRUE} if the data shall be plotted on top of the
#' fitted intensity.
#' @param covariate Character vector of length one which is name of the
#' covariate to be plotted. Must be an external categorical covariate with at
#' most ten different values.
#' @param trans The transformation applied to the color bar of the intensity
#' fit. Specifies the \code{trans} argument of
#' \code{\link[ggplot2]{scale_color_gradient}}.
#' @param sol Solution of the color network plot.
#' @param select Allows the plot for a single model term to be selected for
#' printing. e.g. if just the plot for the second smooth should be printed to
#' the console, set \code{select = 2}.
#' @return Invisibly returns an object of class \code{ggplot} or a list
#' of \code{ggplot} objects.
#' @import ggplot2
#' @export

plot.gn <- function(x, ..., title = list(), size = list(), color = list(),
                    frame = FALSE) {
  a1_x <- a1_y <- a2_x <- a2_y <- NULL
  g <- ggplot(x$lins)
  if (!frame) {
    g <- g + theme_void()
  } else {
    g <- g + theme_bw()
    if (is.null(title$x)) title$x <- paste("x (in", x$unit$plural, ")")
    if (is.null(title$y)) title$y <- paste("y (in", x$unit$plural, ")")
  }
  if (is.null(size$lines)) size$lines <- 1
  if (is.null(color$lines)) color$lines <- "black"
  g <- g + geom_segment(aes(x = a1_x, y = a1_y, xend = a2_x, yend = a2_y),
                        size = size$lines, color = color$lines,
                        lineend = "round", linejoin = "bevel") +
    labs(x = title$x, y = title$y, title = title$plot) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  print(g)
  invisible(g)
}

#' @rdname plot.gn
#' @import ggplot2
#' @export

plot.gnpp <- function(x, ..., title = list(), size = list(), color = list(),
                      shape = 1, frame = FALSE, covariate = NULL) {
  if (length(covariate) > 1){
    stop("Currently, only one covariate can be plotted.")
  }
  y <- a1_x <- a1_y <- a2_x <- a2_y <- NULL
  g <- ggplot(x$network$lins)
  if (!frame) {
    g <- g + theme_void()
  } else {
    g <- g + theme_bw()
    if (is.null(title$x)) title$x <- paste("x (in", x$unit$plural, ")")
    if (is.null(title$y)) title$y <- paste("y (in", x$unit$plural, ")")
  }
  if (is.null(size$lines)) size$lines <- 1
  if (is.null(size$points)) size$points <- 2.5
  if (is.null(color$lines)) color$lines <- "black"
  if (is.null(color$points)) color$points <- "black"
  g <- g + geom_segment(aes(x = a1_x, y = a1_y, xend = a2_x, yend = a2_y),
                        size = size$lines, lineend = "round", linejoin = "mitre") +
    labs(x = title$x, y = title$y, title = title$plot) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  if (is.null(covariate)) {
    g <- g + geom_point(data = x$data, aes(x = x, y = y), shape = shape,
                        color = color$points, size = size$points)
  } else {
    if (is.factor(x$data[[covariate]]) | is.character(x$data[[covariate]])) {
      k <- length(unique(x$data[[covariate]]))
      if (k > 10) stop("Number of different values must be at most ten.")
      g <- g + geom_point(data = x$data, aes(x = x, y = y,
                                             color = !!sym(covariate),
                                             shape = !!sym(covariate)),
                          size = size$points) +
        scale_shape_manual(values = 0:(k-1))
    } else {
      stop("Currently only categorical covariates are supported.")
    }
  }
  print(g)
  invisible(g)
}

#' @rdname plot.gn
#' @import ggplot2
#' @importFrom  grDevices devAskNewPage
#' @importFrom splines splineDesign
#' @export

plot.gnppfit <- function(x, ..., title = list(), size = list(), color = list(),
                         shape = 1, frame = FALSE, data = FALSE,
                         trans = "identity", select = NULL,  sol = 100) {

  stopifnot(inherits(x, "gnppfit"))
  e <- xend <- y <- yend <- intensity <- lower <- upper <- NULL
  G <- as_gn(x)
  g <- df <- vector("list", length(x$smooth) + 1)

  for (m in 1:G$M) {
    lins_m <- filter(G$lins, e == m)
    cs <- c(0, cumsum(lins_m$length))
    dx <- lins_m$a2_x - lins_m$a1_x
    dy <- lins_m$a2_y - lins_m$a1_y
    for (i in 1:length(lins_m$l)) {
      tt <- seq(0, 1, 1/sol)
      xx <- lins_m$a1_x[i] + tt*dx[i]
      yy <- lins_m$a1_y[i] + tt*dy[i]
      zz <- cs[i] + (tt - 1/(2*sol))[-1]*lins_m$length[i]
      df[[1]] <- bind_rows(df[[1]], tibble(l = lins_m$l[i], e = m,
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
    if (is.null(title$x)) title$x <- paste("x (in", x$unit$plural, ")")
    if (is.null(title$y)) title$y <- paste("y (in", x$unit$plural, ")")
  }
  if (is.null(size$lines)) size$lines <- 1
  g[[1]] <- g[[1]] + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = intensity),
                                  size = size$lines, lineend = "round", linejoin = "mitre") +
    labs(color = "Intensity", x = title$x[1], y = title$y[1], title = title$plot[1]) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red", trans = trans)
  if (data) {
    if (is.null(size$points)) size$points <- 2.5
    if (is.null(color$points)) color$points <- "black"
    g[[1]] <- g[[1]] + geom_point(data = X$data, aes(x = x, y = y),
                                  size = size$points, color = color$points,
                                  shape = shape)
  }

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

      g[[i]] <- ggplot(df[[i]]) +
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "grey50") +
        geom_line(aes(x = x, y = y), color = "red") +
        labs(x = title$x[i], y = title$y[i], title = title$plot[i]) +
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
  if (length(g) == 1) g <- g[[1]]
  invisible(g)
}

#' @rdname plot.gn
#' @import ggplot2
#' @export

plot.lppfit <- function(x, ..., title = list(), size = list(), color = list(),
                        shape = 1, frame = FALSE, data = FALSE,
                        trans = "identity", sol = 100) {
  e <- xend <- y <- yend <- intensity <- lower <- upper <- NULL
  G <- as_gn(x)
  f <- as.linfun(x)
  df <- NULL

  for (m in 1:G$M) {
    lins_m <- filter(G$lins, e == m)
    cs <- c(0, cumsum(lins_m$length))
    dx <- lins_m$a2_x - lins_m$a1_x
    dy <- lins_m$a2_y - lins_m$a1_y
    for (i in 1:length(lins_m$l)) {
      tt <- seq(0, 1, 1/sol)
      xx <- lins_m$a1_x[i] + tt*dx[i]
      yy <- lins_m$a1_y[i] + tt*dy[i]
      zz <- cs[i] + (tt - 1/(2*sol))[-1]*lins_m$length[i]
      df <- bind_rows(df, tibble(l = lins_m$l[i], e = m,
                                 x = utils::head(xx, -1), xend = xx[-1],
                                 y = utils::head(yy, -1), yend = yy[-1],
                                 z = zz))
    }
  }
  df$intensity <- f(df$x, df$y)
  if (!frame) {
    g <- ggplot(df) + theme_void()
  } else {
    g <- ggplot(df) + theme_bw()
    if (is.null(title$x)) title$x <- paste("x (in", x$unit$plural, ")")
    if (is.null(title$y)) title$y <- paste("y (in", x$unit$plural, ")")
  }
  if (is.null(size$lines)) size$lines <- 1
  g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = intensity),
                    size = size$lines, lineend = "round", linejoin = "mitre") +
    labs(color = "Intensity", x = title$x, y = title$y, title = title$title) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red", trans = trans)
  if (data) {
    if (is.null(size$points)) size$points <- 2.5
    if (is.null(color$points)) color$points <- "black"
    g <- g + geom_point(data = x$data, aes(x = x, y = y),
                                  size = size$points, color = color$points,
                                  shape = shape)
  }
  print(g)
  invisible(g)
}
