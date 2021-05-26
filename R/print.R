#' Print Method for Geometric Networks
#'
#' Prints basic information of a geometric network related object.
#'
#' @param x A geometric network (object of class \code{gn}), a
#' point pattern on a geometric network (object of class \code{gnpp}) or a
#' fitted point process on a geometric network (object of class \code{gnppfit}).
#' @param ... Other arguments.
#' @export

print.gn <- function(x, ...){
  stopifnot(inherits(x, "gn"))
  cat(paste("Geometric network in", x$q, "dimensions with",
            x$W, "vertices and", x$M, "line segments.\n"))
  cat(paste("Length of the network:", round(sum(x$d), 3), "units.\n"))
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
  cat(paste("Length of the network:", round(sum(x$d), 3), "units.\n"))
  if (ncol(x$data) > 6) {
    n <- ncol(x$data) - 6
    if (n == 1){
      cat("Point hattern has one covariate:\n")
    } else {
      cat(paste0("Point pattern has ", n, " covariates:\n"))
    }

    for (i in 1:n) {
      covariate <- x$data[[i + 6]]
      cl <- as.character(class(covariate))
      cat(paste0(i, ") '", names(x$data[i + 6]), "': Covariate of class ",
                 cl, " "))
      if (cl == "numeric") {
        cat(paste0("with ",  length(unique(covariate)), " unique values.\n"))
      }
      if (cl == "factor") {
        cat(paste0("with ", nlevels(covariate), " levels.\n"))
      }
      summary_cov <- t(as.matrix(summary(x$data %>% pull(7))))
      rownames(summary_cov) <- ""
      printCoefmat(summary_cov)
      cat("\n")
    }
  } else {
    cat("Point pattern has no covariates.")
  }
  invisible(NULL)
}

#' @rdname print.gn
#' @export

print.gnppfit <- function(x, ...){
  stopifnot(inherits(x, "gnppfit"))
  cat(paste("Geometric network in", x$network$q, "dimensions"))
  invisible(NULL)
}

#' Print Method for gnppfit Summary
#'
#' Prints basic information on the model fit and the underlying geometric
#' network.
#'
#' @param x A model fit on a geometric network (object of class \code{gnpp}).
#' @param ... further arguments passed to print
#' @importFrom stats printCoefmat
#' @export

print.summary.gnppfit <- function(x, ...){
  cat("Poisson model on a geometric network fitted with maximum likelihood")
  cat("\n\nFormula:\n")
  print(x$formula)
  cat("\nParametric coefficients:\n")
  printCoefmat(x$tab[, c(1:2, 6:7)], #digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  cat(paste0("\nNumber of Fellner-Schall-iterations: ", x$it_rho))
}
