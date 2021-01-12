## code to prepare `test` dataset goes here
library(dplyr)
small_geonet <- list(
  vertices = tibble(x = c(0, 1, 2, 2, 2, 3, 1), y = c(0, 1, 3, 2, 0, 3, 2), id = 1:7),
  curves = list(
    tibble(from = 1, to = 2, m = 1),
    tibble(from = 2, to = 4, m = 2),
    tibble(from = c(2, 7), to = c(7, 3), m = 3),
    tibble(from = 3, to = 4, m = 4),
    tibble(from = 4, to = 5, m = 5),
    tibble(from = 4, to = 6, m = 6)
  )
)

class(small_geonet) <- "geonet"
usethis::use_data(small_geonet, overwrite = TRUE)
