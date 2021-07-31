## code to prepare `small_geonet` dataset goes here
x <- list(
  vertices = dplyr::tibble(
    a = 1:7,
    v = c(1, 2, 3, 4, 5, 6, NA),
    x = c(0, 1, 2, 2, 2, 3, 1),
    y = c(0, 1, 3, 2, 0, 3, 2)
  ),
  lins = dplyr::tibble(
    l = 1:7,
    e = c(1, 2, 3, 3, 4, 5, 6),
    a1 = c(1, 2, 2, 7, 3, 4, 4),
    a2 = c(2, 4, 7, 3, 4, 5, 6)
  ),
  adjacency = matrix(0, 6, 6),
  incidence = matrix(0, 6, 6),
  d = rep(0, 6),
  q = 2,
  W = 6,
  M = 6
)
x$lins$a1_x <- x$vertices$x[x$lins$a1]
x$lins$a1_y <- x$vertices$y[x$lins$a1]
x$lins$a2_x <- x$vertices$x[x$lins$a2]
x$lins$a2_y <- x$vertices$y[x$lins$a2]
x$lins$length <- sqrt((x$lins$a1_x - x$lins$a2_x)^2 +
  (x$lins$a1_y - x$lins$a2_y)^2)
x$d <- x$lins %>%
  dplyr::group_by(e) %>%
  dplyr::summarize(length = sum(length), .groups = "drop") %>%
  dplyr::pull(length)
for (m in 1:x$M) {
  dat <- dplyr::filter(x$lins, e == m)
  x$adjacency[dat$a1[1], dat$a2[nrow(dat)]] <- 1
  cs <- cumsum(dplyr::filter(x$lins, e == m) %>% dplyr::pull(length))
  x$lins$frac1[which(x$lins$e == m)] <- c(0, cs[-length(cs)]/cs[length(cs)])
  x$lins$frac2[which(x$lins$e == m)] <- dat$length/sum(dat$length)
}
x$incidence <- incidence(x$vertices, x$lins)
x$unit <- unitname(c("unit", "units"))

small_gn <- x
class(small_gn) <- "gn"
usethis::use_data(small_gn, overwrite = TRUE)
