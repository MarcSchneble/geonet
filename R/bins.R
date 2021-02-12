#' Internal Covariates
#'
#' \code{internal} computes the values of internal covariates in every
#' bin of the network.
#'
#' @param vars The name of the covariates which should go into the model
#' as linear internal covariates.
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param P Penalized spline representation of the geometric network \code{X}.
#' @return A data frame with the number of rows as the number of bins on the
#' geometric network (\code{sum(P$bins$N)}) and the number of columns equal to
#' the length of \code{vars}.
#' @export

internal <- function(vars, X, P){
  out <- list()
  if (is.element("x", vars)) {
    x <- rep(NA, sum(P$bins$N))
    ind <- 1
    for (m in 1:X$network$M) {
      lins_m <- filter(X$network$lins, e == m)
      tp <- P$bins$z[[m]]/X$network$d[m]
      dx <- lins_m$v2_x - lins_m$v1_x
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        x[ind:(ind + length(ind_id) - 1)] <- lins_m$v1_x[i] + tp[ind_id]*dx[i]
        ind <- ind + length(ind_id)
      }
    }
    out$x <- x/1000
  }
  if (is.element("y", vars)) {
    y <- rep(NA, sum(P$bins$N))
    ind <- 1
    for (m in 1:X$network$M) {
      lins_m <- filter(X$network$lins, e == m)
      tp <- P$bins$z[[m]]/X$network$d[m]
      dy <- lins_m$v2_y - lins_m$v1_y
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        y[ind:(ind + length(ind_id) - 1)] <- lins_m$v1_y[i] + tp[ind_id]*dy[i]
        ind <- ind + length(ind_id)
      }
    }
    out$y <- y/1000
  }
  out
}

#' Bin Point Pattern on a Geometric Network
#'
#' \code{binData} bins the data on the supplied point pattern according to
#' all possible combination of covariates.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param P as
#' @param lins as
#' @param smooths as
#' @return The binned data.
#' @import dplyr
#' @importFrom stats setNames
#' @export

binData <- function(X, P, vars, intern){
  # get covariates from every point on the network (if applicable)
  vars_external <- setdiff(vars, c("G", "x", "y"))
  if (all(vars_external %in% colnames(X$data)[-(1:5)])) {
    covariates <- as_tibble(X$data) %>% select(all_of(vars_external))
  } else {
    stop("At least one covariate was not found in data!")
  }
  if (ncol(covariates) == 0) covariates <- covariates %>%mutate(a = 1)

  # get all combinations of covariates and calculate the number of rows of the data matrix
  vars_comb <- covariates %>% distinct() %>% expand.grid() %>%
    distinct() %>% as_tibble()

  # sort data frame
  for (a in 1:length(vars_comb)) {
    vars_comb <- arrange(vars_comb, !!sym(names(vars_comb[a])))
  }

  # initializing
  data <- tibble(id = 1:sum(P$bins$N), count = NA, h = NA)
  if (length(intern) > 0) {
    data <- bind_cols(data, as_tibble(internal(vars[intern], X, P))) %>%
      slice(rep(1:n(), nrow(vars_comb)))
  }
  data <- bind_cols(data, vars_comb %>% slice(rep(1:n(), each = sum(P$bins$N))))

  ind <- 1
  for (j in 1:nrow(vars_comb)) {
    ind_comb <- which(do.call(paste, covariates) == do.call(paste, vars_comb[j, ]))
    data_sub <- X$data[ind_comb, ]
    for (m in 1:X$network$M) {
      # positions of data on curve e
      ind_e <- which(data_sub$e == m)
      y_e <- sort(as.numeric(data_sub[ind_e, ]$tp))*X$network$d[m]

      # bin data
      y_b <- rep(0, length(P$bins$z[[m]]))
      for (k in 1:length(P$bins$z[[m]])) {
        y_b[k] <- length(which(y_e < P$bins$b[[m]][k+1] & y_e > P$bins$b[[m]][k]))
      }

      # stack into one vector
      data$count[ind:(ind + length(y_b) - 1)] <- y_b
      data$h[ind:(ind + length(y_b) - 1)] <- P$bins$h[m]
      ind <- ind + length(y_b)
    }
    # add bin id for every row
    data$id[((j-1)*sum(P$bins$N) + 1):(j*sum(P$bins$N))] <- 1:sum(P$bins$N)
  }
  data$offset <- 1
  data
}
