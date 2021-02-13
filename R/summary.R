#' Intensity Estimation on Geometric Networks with Penalized Splines
#'
#' \code{intensityPspline} estimates the intensity of a point pattern on a
#' geometric network.
#'
#' @param object Point pattern on a geometric network (object of class \code{gnpp})
#' @param ... asdad
#' @return A fitted geometric network (object of class \code{gnppfit}).
#' @import dplyr
#' @export
#'

summary.gnppfit <- function(object, ...) {

  if (length(object$ind[["lins"]]) > 0) {
    tab <- tibble(name = names(object$ind[["lins"]]),
                  estimate = NA, se = NA, rr = NA, rr.lower = NA, rr.upper = NA)
    for (i in 1:length(object$ind[["lins"]])) {
      tab$estimate[i] <- round(object$coefficients[object$ind[["lins"]][i]], 3)
      tab$se[i] <- round(sqrt(object$V[object$ind[["lins"]][i], object$ind[["lins"]][i]]), 3)
      tab$rr[i] <- round(exp(tab$estimate[i]), 2)
      tab$rr.lower[i] <- round(exp(tab$estimate[i] - 1.96*tab$se[i]), 2)
      tab$rr.upper[i] <- round(exp(tab$estimate[i] + 1.96*tab$se[i]), 2)
    }
    tab
  } else {
    "No linear effects to display."
  }
}
