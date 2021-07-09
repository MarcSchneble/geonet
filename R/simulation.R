#' Simulate Uniform Points on a Geometric Network
#'
#' \code{runifgn} simulates points on a geometric network according
#' to a uniform density.
#'
#' @param n Number of random points to generate. A nonnegative integer.
#' @param G A geometric network (object of class \code{gn}).
#' @return An object of class \code{gnpp}.
#' @import dplyr
#' @importFrom stats runif
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export


runifgn <- function(n, G){
  e <- seg <- tp <- frac1 <- frac2 <- x <- y <- NULL
  lins <- G$lins %>% arrange(id)
  cumlen <- c(0, cumsum(lins$length))
  dx <- lins$v2_x - lins$v1_x
  dy <- lins$v2_y - lins$v1_y
  uu <- runif(n, min = 0, max = sum(G$d))
  ii <- findInterval(uu, cumlen, rightmost.closed = TRUE, all.inside = TRUE)
  tt <- (uu - cumlen[ii])/lins$length[ii]
  xx <- lins$v1_x[ii] + tt*dx[ii]
  yy <- lins$v1_y[ii] + tt*dy[ii]
  data <- tibble(id = ii, tp = tt)
  data <- left_join(data, lins, by = "id") %>%
    mutate(x = xx, y = yy, tp = frac1 + tp*frac2) %>%
    select(id, e, tp, x, y) %>%
    arrange(e, tp)
  out <- list(data = data, network = G)
  class(out) <- "gnpp"
  out
}

rgnpp <- function(n, fit){
  int <- vector("list", fit$network$M)
  for (m in 1:fit$network$M) {
    for (k in 1:(fit$knots$J[m]+1)) {
      int[[m]][k] <- integrate(network_intensity_edge,
                              lower = fit$knots$tau[[m]][k],
                              upper = fit$knots$tau[[m]][k+1],
                              m = m, fit1 = fit)$value
    }
  }
  intens_total <- sum(unlist(int))
  prob_edges <- unlist(lapply(int, function(x) sum(x)))/intens_total
  mm <- table(sample(1:fit$network$M, size = n, replace = TRUE, prob = prob_edges))
  i <- 0
  for (m in as.numeric(names(mm))) {
    i <- i + 1

    ind_v1 <- which(fit$network$incidence[, m] == -1)
    ind_v2 <- which(fit$network$incidence[, m] == 1)
    gamma <- as.numeric(fit$coefficients[fit$ind[[1]]])
    gamma_m <- c(gamma[sum(fit$knots$J) + ind_v1],
                  gamma[((cumsum(fit$knots$J)-fit$knots$J)[m]+1):cumsum(fit$knots$J)[m]],
                  gamma[sum(fit$knots$J) + ind_v2])

    probs <- int[[m]]/sum(int[[m]])
    kk <- sample(1:(fit$knots$J[m]+1), size = mm[i], replace = TRUE, prob = probs)

    uu <- runif(mm[i], min = 0, max = 1)

    z <- log(uu*int[[m]][kk]*(gamma_m[kk+1] - gamma_m[kk])/
               (exp(gamma_m[kk])*fit$knots$delta[m]) + 1)*
      fit$knots$delta[m]/(gamma_m[kk+1] - gamma_m[kk])

  }
}
