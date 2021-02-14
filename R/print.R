#' Print Method for Geometric Networks
#'
#' Prints basic information on the geometric network
#'
#' @param x A geometric network (object of class \code{gn}) or a
#' point pattern on a geometric network (object of class \code{gnpp}).
#' @param ... further arguments passed to print
#' @export

print.gn <- function(x, ...){
  stopifnot(inherits(x, "gn"))
  cat(paste("Geometric network in", x$q, "dimensions with",
            x$W, "vertices and", x$M, "line segments."))
  invisible(NULL)
}

#' @rdname print.gn
#' @export

print.gnpp <- function(x, ...){
  stopifnot(inherits(x, "gnpp"))
  cat(paste0("Point pattern of size ", nrow(x$data),
            " on a geometric network in ",
            x$network$q, " dimensions.\n",
            "Network has ", x$network$W, " vertices and ",
            x$network$M, " curve segments.\n"))
  if (ncol(x$data) > 5) {
    n <- ncol(x$data) - 5
    if (n == 1){
      cat("Point hattern has one covariate:\n")
    } else {
      cat(paste0("Point pattern has ", n, " covariates:\n"))
    }

    for (i in 1:n) {
      covariate <- x$data[[i + 5]]
      cl <- as.character(class(covariate))
      cat(paste0(i, ") '", names(x$data[i + 5]), "': Covariate of class ",
                 cl, " "))
      if (cl == "numeric") {
        cat(paste0("with minimum value ", min(covariate),
                   ", maximum value ", max(covariate),
                   " and ", length(unique(covariate)), " unique values."))
      }
      if (cl == "factor") {
        cat(paste0("with ", nlevels(covariate), " levels."))
      }
    }
  } else {
    cat("Point pattern has no covariates.")
  }
  invisible(NULL)
}

#' Print Method for a fitted Geometric Network Model
#'
#' Prints basic information on the model fit and the underlying geometric
#' network.
#'
#' @param x A model fit on a geometric network (object of class \code{gnpp}).
#' @param ... further arguments passed to print
#' @export
#'
print.gnppfit <- function(x, ...){
  stopifnot(inherits(x, "gnppfit"))
  cat(paste("Geometric network in", x$network$q, "dimensions"))
  invisible(NULL)
}
