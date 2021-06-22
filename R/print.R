#' Print Method for Geometric Networks
#'
#' Prints basic information of a geometric network related object.
#'
#' @param x A geometric network (object of class \code{gn}), a
#' point pattern on a geometric network (object of class \code{gnpp}) or a
#' fitted point process on a geometric network (object of class \code{gnppfit}).
#' @param ... Other arguments.
#' @return Invisibly returns the NULL object.
#' @export

print.gn <- function(x, ...){
  stopifnot(inherits(x, "gn"))
  cat(paste("Geometric network in", x$q, "dimensions with",
            x$W, "vertices and", x$M, "curve segments.\n"))
  invisible(NULL)
}

#' @rdname print.gn
#' @importFrom stats printCoefmat
#' @importFrom dplyr pull
#' @export

print.gnpp <- function(x, ...){
  stopifnot(inherits(x, "gnpp"))
  cat(paste0("Point pattern of size ", nrow(x$data),
             " on a geometric network in ",
             x$network$q, " dimensions.\n",
             "Network has ", x$network$W, " vertices and ",
             x$network$M, " curve segments.\n"))
  cat("Point pattern has ")
  if (ncol(x$network$lins) == 11) {
    cat("no internal covariates and ")
  } else {
    if (ncol(x$network$lins) == 12) {
      cat("1 internal covariate and ")
    } else {
      cat(paste(ncol(x$network$lins) - 11), "internal covariates and ")
    }
  }
  if (ncol(x$data) == 6) {
    cat("no external covariates.")
  } else if (ncol(x$data) == 7) {
    cat("1 external covariate.")
  } else {
    cat(paste(ncol(x$data) - 6, "external covariates."))
  }
  invisible(NULL)
}

#' @rdname print.gn
#' @export

print.gnppfit <- function(x, ...){
  stopifnot(inherits(x, "gnppfit"))
  cat(paste0("Fitted point prozess on a geometric network in ",
             x$network$q, " dimensions.\n",
             "Network has ", x$network$W, " vertices and ",
             x$network$M, " curve segments.\n"))
  invisible(NULL)
}

#' Print Method for Summaries
#'
#' Prints basic information of a geonet summary object.
#'
#' @param x A geonet summary object of class \code{summary.gn},
#' code{summary.gnpp} or \code{summary.gnppfit}.
#' @param ... Other arguments.
#' @return Invisibly returns the NULL object.
#' @export

print.summary.gn <- function(x, ...) {
  stopifnot(inherits(x, "summary.gn"))
  cat(paste("Geometric network in", x$q, "dimensions with",
            x$W, "vertices and", x$M, "curve segments.\n"))
  cat(paste("Total length of the network:", round(x$total_length, 3),
            x$unit$plural, "\n"))
  cat(paste("Minimum segment length:", round(x$range_length[1], 3),
            x$unit$plural, "\n"))
  cat(paste("Maximum segment length:", round(x$range_length[2], 3),
            x$unit$plural, "\n"))
  cat(paste("Distribution of vertex degrees:"))
  print(x$degrees)
}

#' @rdname print.summary.gn
#' @export

print.summary.gnpp <- function(x, ...) {
  stopifnot(inherits(x, "summary.gnpp"))
  cat(paste("Point pattern on a geometric network in", x$q, "dimensions with",
            x$W, "vertices and", x$M, "curve segments.\n"))
  cat(paste("Total length of the network:", round(x$total_length, 3),
            x$unit$plural, "\n"))
  if (length(x$covariates$internal) > 0) {
    cat(paste("Number of network internal covariates:",
              length(x$covariates$internal), "\n"))
    for (i in 1:length(x$covariates$internal)) {
      cat(paste(i, ") ", x$covariates$internal[[i]]$class, " variable \"",
                names(x$covariates$internal)[i], "\":\n", sep = ""))
      print(x$covariates$internal[[i]]$summary)
    }
  } else {
    cat("Network has no internal covariates\n")
  }
  cat(paste("Number of points:", x$n_points, "\n"))
  cat(paste("Average intensity:", round(x$average_intensity, 5), "per",
            x$unit$singular, "\n"))
  if (length(x$covariates$internal) > 0) {
    cat(paste("Number of external covariates:",
              length(x$covariates$external), "\n"))
    for (j in (i+1):(length(x$covariates$external)+i)) {
      cat(paste(j, ") ", x$covariates$external[[j-i]]$class, " variable \"",
                names(x$covariates$external)[j-i], "\":\n", sep = ""))
      print(x$covariates$external[[j-i]]$summary)
    }
  } else {
    cat("Point pattern has no external covariates\n")
  }
}

#' @rdname print.summary.gn
#' @importFrom stats printCoefmat
#' @export

print.summary.gnppfit <- function(x, ...){
  cat("Intensity estimation on a geometric network.\n")
  cat("Family: Poisson\n")
  cat("Link Function: log\n")
  cat("\nFormula: ")
  print(x$formula, showEnv = FALSE)
  if (!is.null(x$tab)) {
    cat("\nPparametric coefficients:\n")
    printCoefmat(x$tab[, c(1:2, 6:7)], #digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  } else {
    cat("\nModel has no parametric coefficients.\n")
  }
  cat("\nApproximate significance of smooth terms:\n")
  print(x$edf)
  cat(paste0("\nNumber of Fellner-Schall-iterations: ", x$it_rho))
  invisible(NULL)
}
