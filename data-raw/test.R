## code to prepare `small_geonet` dataset goes here
library(dplyr)
x <- list(
  vertices = tibble(
    id = 1:7,
    x = c(0, 1, 2, 2, 2, 3, 1),
    y = c(0, 1, 3, 2, 0, 3, 2)
  ),
  lins = tibble(
    id = 1:7,
    m = c(1, 2, 3, 3, 4, 5, 6),
    v1 = c(1, 2, 2, 7, 3, 4, 4),
    v2 = c(2, 4, 7, 3, 4, 5, 6)
  ),
  d = rep(0, 6),
  A_v = matrix(0, 6, 6),
  q = 2,
  W = 6,
  M = 6
)
x$lins$v1_x <- x$vertices$x[x$lins$v1]
x$lins$v1_y <- x$vertices$y[x$lins$v1]
x$lins$v2_x <- x$vertices$x[x$lins$v2]
x$lins$v2_y <- x$vertices$y[x$lins$v2]
x$lins$length <- sqrt((x$lins$v1_x - x$lins$v2_x)^2 +
  (x$lins$v1_y - x$lins$v2_y)^2)
x$d <- x$lins %>% group_by(m) %>% summarize(length = sum(length)) %>% pull(length)
for (i in 1:x$M) {
  sub <- filter(x$lins, m == i)
  x$A_v[sub$v1[1], sub$v2[nrow(sub)]] <- 1
}
x$A_v <- pmin(x$A_v + t(x$A_v), 1)

small_gn <- x
class(small_gn) <- "gn"
usethis::use_data(small_gn, overwrite = TRUE)
