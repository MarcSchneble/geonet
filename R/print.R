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

#' Print Method for Geometric Networks Data
#'
#' This is the print method for object of the class gnd.
#'
#' @param x an object of class gns
#' @param ... further arguments passed to print
#' @export
print.gnpp <- function(x, ...){
  stopifnot(inherits(x, "gnpp"))
  cat(paste("Geometric network in", x$network$q, "dimensions"))
  invisible(NULL)
}

#' Print Method for Geometric Networks Splines
#'
#' This is the print method for object of the class gns.
#'
#' @param x an object of class gns
#' @param ... further arguments passed to print
#' @export
print.gnppfit <- function(x, ...){
  stopifnot(inherits(x, "gnppfit"))
  cat(paste("Geometric network in", x$network$q, "dimensions"))
  invisible(NULL)
}
