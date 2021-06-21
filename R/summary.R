#' Summary for a (point pattern on a) geometric network
#'
#' Takes an object of class \code{gn} or \code{gnpp} and computes a summary
#' from it.
#'
#' @param object A geometric network (object of class \code{gn} or a point
#' pattern on a geometric network (object of class \code{gnpp}).
#' @param ... Other arguments.
#' @return A summary of a (point pattern on a) geometric network. This is an
#' object of class \code{summary.gn} or \code{summary.gnpp}, respectively.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

summary.gn <- function(object, ...) {
  stopifnot(inherits(object, "gn"))
  deg <- colSums(object$adjacency)
  out <- list(
    M = object$M,
    W = object$W,
    q = object$q,
    unit = object$unit,
    total_length = sum(object$d),
    range_length = range(object$d),
    degrees = table(colSums(object$adjacency))
  )
  class(out) <- c("summary.gn", class(out))
  out
}

#' @rdname summary.gn
#' @export

summary.gnpp <- function(object, ...) {
  stopifnot(inherits(object, "gnpp"))
  G <- as_gn(object)
  out <- summary(G)
  out$n_points <- nrow(object$data)
  out$average_intensity <- out$n_points/out$total_length
  out$covariates <- list(internal = NULL, external = NULL)
  if (ncol(G$lins) > 11) {
    nam <- colnames(G$lins)[-(1:11)]
    out$covariates$internal <- vector("list", ncol(G$lins) - 11)
    names(out$covariates$internal) <- nam
    for (i in 1:length(nam)) {
      out$covariates$internal[[i]] <- list(
        class = class(G$lins[[paste(nam[i])]]),
        summary = summary(G$lins[[paste(nam[i])]])
      )
    }
  }
  if (ncol(object$data) > 6) {
    nam <- colnames(object$data)[-(1:6)]
    out$covariates$external <- vector("list", ncol(object$data) - 6)
    names(out$covariates$external) <- nam
    for (i in 1:length(nam)) {
      out$covariates$external[[i]] <- list(
        class = class(object$data[[paste(nam[i])]]),
        summary = summary(object$data[[paste(nam[i])]]))
    }
  }
  class(out) <- c("summary.gnpp", class(out))
  out
}

#' Summary for fitted point process on a geometric network
#'
#' Takes a fitted \code{gnppfit} object produced by \code{intensity_spline}
#' and computes a summary from it.
#'
#' @param object A fitted point process on a geometric network.
#' @param ... Other arguments.
#' @return A summary of a fitted point process on a geometric network. This is
#' an object of class \code{summary.gnppfit}.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @importFrom stats pnorm
#' @export

summary.gnppfit <- function(object, ...) {

  x <- object
  ans <- list()
  if (length(x$ind[["lins"]]) > 0) {
    ans$tab <- data.frame(estimate = numeric(length(x$ind[["lins"]])),
                          se = NA, rr = NA, rr.lower = NA, rr.upper = NA,
                          z = NA, p = NA)
    for (i in 1:length(x$ind[["lins"]])) {
      ans$tab$estimate[i] <- x$coefficients[x$ind[["lins"]][i]]
      ans$tab$se[i] <- sqrt(x$V[x$ind[["lins"]][i], x$ind[["lins"]][i]])
      ans$tab$rr[i] <- exp(ans$tab$estimate[i])
      ans$tab$rr.lower[i] <- exp(ans$tab$estimate[i] - 1.96*ans$tab$se[i])
      ans$tab$rr.upper[i] <- exp(ans$tab$estimate[i] + 1.96*ans$tab$se[i])
      ans$tab$z[i] <- ans$tab$estimate[i]/ans$tab$se[i]
      ans$tab$p[i] <- (1 - pnorm(abs(ans$tab$z[i])))*2
    }
    colnames(ans$tab) <- c("Estimate", "Std. Error", "rr", "rr.lower",
                           "rr.upper", "z value", "Pr(>|z|)")
    rownames(ans$tab) <- names(x$ind[["lins"]])
  }
  edf <- rep(NA, length(x$ind))
  for (k in 1:length(x$ind)) {
    edf[k] <- sum(x$edf[x$ind[[k]]])
  }
  if (is.element("lins", names(x$ind))) {
    edf <- edf[1:(length(edf) - 1)]
  }
  ans$edf <- edf
  ans$formula <- x$formula
  ans$it_rho <- x$it_rho
  class(ans) <- "summary.gnppfit"
  ans
}
