#' Methods for Geometric Networks
#'
#' The function as.gns converts an object of class gn to an object of
#' class gns
#'
#' @param x an object of class gn or an object of class linnet (or an object that can
#' be converted to an instance of class linnet)
#' @return an object of class gn
ends <- function(x){
  ends <- x$lins %>% filter(m == 0)
  for (i in 1:x$M) {
    sub <- x$lins %>% filter(m == i)
    ends <- bind_rows(ends, sub %>% slice(1))
    ends$to[i] <- sub$to[nrow(sub)]
    ends$length[i] <- sum(sub$length)
  }
  ends
}
