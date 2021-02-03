getIncidence <- function(x){
  incidence <- matrix(0, x$W, x$M)
  for (m in 1:x$M) {
    dat <- dplyr::filter(x$lins, e == m)
    incidence[x$vertices$v[match(dat$v1[1], x$vertices$id)], m] <- -1
    incidence[x$vertices$v[match(tail(dat$v2, 1), x$vertices$id)], m] <- 1
  }
  incidence
}
