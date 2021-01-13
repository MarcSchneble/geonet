## code to prepare `small_geonet` dataset goes here
library(dplyr)
x <- list(
  q = 2,
  vertices = tibble(
    id = 1:7,
    x = c(0, 1, 2, 2, 2, 3, 1),
    y = c(0, 1, 3, 2, 0, 3, 2)
  ),
  lins = tibble(
    m = c(1, 2, 3, 3, 4, 5, 6),
    from = c(1, 2, 2, 7, 3, 4, 4),
    to = c(2, 4, 7, 3, 4, 5, 6),
    length = NA
  ),
  curves = list(),
  d = rep(0, 6),
  d_L = 0,
  M = 6,
  W = 6,
  A_v = matrix(0, 6, 6)
)
x$lins$length <- sqrt((x$vertices$x[x$lins$from] - x$vertices$x[x$lins$to])^2 +
  (x$vertices$y[x$lins$from] - x$vertices$y[x$lins$to])^2)
x$d <- x$lins %>% group_by(m) %>% summarize(length = sum(length)) %>% pull(length)
x$d_L <- sum(x$d)
for (i in 1:x$M) {
  sub <- filter(x$lins, m == i)
  x$A_v[sub$from[1], sub$to[nrow(sub)]] <- 1
}
x$A_v <- pmin(x$A_v + t(x$A_v), 1)

small_gn <- x
class(small_gn) <- "gn"
usethis::use_data(small_gn, overwrite = TRUE)
