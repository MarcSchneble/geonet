#' Simulate Uniform Points on a Geometric Network
#'
#' \code{runifgn} simulates points on a geometric network according
#' to a uniform density.
#'
#' @param n Number of random points to generate. A nonnegative integer.
#' @param G A geometric network (object of class \code{gn}).
#' @return A point pattern on a geometric network, an object of class
#' \code{gnpp}.
#' @import dplyr
#' @importFrom stats runif
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export


runifgn <- function(n, G){
  l <- tp_e <- e <- seg <- tp_l <- frac1 <- frac2 <- x <- y <- NULL
  lins <- G$lins %>% arrange(l)
  cumlen <- c(0, cumsum(lins$length))
  dx <- lins$a2_x - lins$a1_x
  dy <- lins$a2_y - lins$a1_y
  uu <- runif(n, min = 0, max = sum(G$d))
  ii <- findInterval(uu, cumlen, rightmost.closed = TRUE, all.inside = TRUE)
  tt <- (uu - cumlen[ii])/lins$length[ii]
  xx <- lins$a1_x[ii] + tt*dx[ii]
  yy <- lins$a1_y[ii] + tt*dy[ii]
  data <- tibble(l = ii, tp_l = tt)
  data <- left_join(data, lins, by = "l") %>%
    mutate(x = xx, y = yy, tp_e = frac1 + tp_l*frac2) %>%
    select(l, tp_l, e, tp_e, x, y) %>%
    arrange(e, tp_l)
  out <- list(data = data, network = G)
  class(out) <- "gnpp"
  out
}

#' Random Points on a Geometric Network
#'
#' @param n Number of random points. A nonnegative integer.
#' @param fit A fitted point process on a geometric network (object of class
#' \code{gnpp}).
#'
#' @return A point pattern on a geometric network, an object of class
#' \code{gnpp}.
#' @importFrom stats integrate
#' @import dplyr
#' @export
#' @examples
#' X <- runifgn(n = 50, G = small_gn)
#' fit <- intensity_pspline(X)
#' X2 <- rgnpp(n = 50, fit = fit)
#' plot(X2)

rgnpp <- function(n, fit){
  e <- tp_l <- NULL
  int <- vector("list", fit$network$M)
  for (m in 1:fit$network$M) {
    for (k in 1:(fit$knots$J[m]+1)) {
      int[[m]][k] <- integrate(network_intensity,
                              lower = fit$knots$tau[[m]][k],
                              upper = fit$knots$tau[[m]][k+1],
                              m = m, fit1 = fit)$value
    }
  }
  intens_total <- sum(unlist(int))
  prob_edges <- unlist(lapply(int, function(x) sum(x)))/intens_total
  mm <- table(sample(1:fit$network$M, size = n, replace = TRUE, prob = prob_edges))
  data <- tibble(l = 1:n, tp_l = numeric(1), e = integer(1), tp_e = numeric(1), x = numeric(1), y = numeric(1))
  ind <- 1
  i <- 0
  for (m in as.numeric(names(mm))) {
    i <- i + 1
    # get the coefficients related to the m-th curve
    ind_v1 <- which(fit$network$incidence[, m] == -1)
    ind_v2 <- which(fit$network$incidence[, m] == 1)
    gamma <- as.numeric(fit$coefficients[fit$ind[[1]]])
    gamma_m <- c(gamma[sum(fit$knots$J) + ind_v1],
                  gamma[((cumsum(fit$knots$J)-fit$knots$J)[m]+1):cumsum(fit$knots$J)[m]],
                  gamma[sum(fit$knots$J) + ind_v2])

    probs <- int[[m]]/sum(int[[m]])
    kk <- sample(1:(fit$knots$J[m]+1), size = mm[i], replace = TRUE, prob = probs)

    uu <- runif(mm[i], min = 0, max = 1)

    # check whether results is between 0 and 1
    z <- log(uu*int[[m]][kk]*(gamma_m[kk+1] - gamma_m[kk])/
               (exp(kk*gamma_m[kk] + gamma_m[kk+1]*(1-kk))*fit$knots$delta[m]) +
               exp((gamma_m[kk+1]-gamma_m[kk])/fit$knots$delta[m]*fit$knots$tau[[m]][kk]))*
      fit$knots$delta[m]/(gamma_m[kk+1] - gamma_m[kk])

    coordinates <- network_location(fit$network, m = m, z = z)
    data$l[ind:(ind+mm[i]-1)] <- coordinates$l
    data$tp_l[ind:(ind+mm[i]-1)] <- coordinates$tp_l
    data$e[ind:(ind+mm[i]-1)] <- m
    data$tp_e[ind:(ind+mm[i]-1)] <- z/fit$network$d[m]
    data$x[ind:(ind+mm[i]-1)] <- coordinates$x
    data$y[ind:(ind+mm[i]-1)] <- coordinates$y
    ind <- ind + as.numeric(mm[i])
    if (min(data$tp_e) < 0 | max(data$tp_e) > 1) {
      warning("Data were simulated outside of the range of the network.")
    }
  }
  data <- data %>% arrange(e, tp_l)
  out <- list(data = data, network = fit$network)
  class(out) <- "gnpp"
  out
}
