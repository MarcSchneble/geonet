#' Summary for fitted point process on a geometric network
#'
#' Takes a fitted \code{gnppfit} object produced by \code{intensity_spline} and
#' computes a summary from it.
#'
#' @param object Fitted point process on a geometric network (object of class
#' \code{gnppfit}).
#' @param ... Other arguments.
#' @return A summary of a fitted point process on a geometric network (object of
#' class \code{summary.gnppfit}).
#' @importFrom stats pnorm
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
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
