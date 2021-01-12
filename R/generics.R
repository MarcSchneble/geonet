#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
print.geonet <- function(object, ...){
  stopifnot(inherits(object, "geonet"))
  cat(paste("Geometric network in", 2, "dimensions"))
  invisible(NULL)
}

#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
plot.geonet <- function(object, ..., x = "x", y = "y", title = NULL){
  dat <- object$curves %>%
    mutate(x.from = object$vertices$x[from],
           y.from = object$vertices$y[from],
           x.to = object$vertices$x[to],
           y.to = object$vertices$y[to])
  g <- ggplot(dat) +
    geom_segment(aes(x = x.from, y = y.from, xend = x.to, yend = y.to)) +
    labs(x = x, y = x, title = title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  print(g)
  invisible(g)
}
