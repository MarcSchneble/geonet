#' Penalized Splines on a Geometric Network
#'
#' \code{pspline} constructs linear B-splines on a geometric network and
#' subdivides each curve segment into severval bins.
#'
#' @param G A geometric network (object of class \code{gn}) or,
#' alternatively, a point pattern on a geometric network
#' (object of class \code{gnpp}).
#' @param delta The knot distance delta.
#' @param h The bin width h.
#' @return Linear B-Splines and bins on the geometric network.
#' @export

pspline <- function(G, delta = NULL, h = NULL){
  if (inherits(G, "gnpp")) G <- as_gn(G)
  if (!inherits(G, "gn")) stop("G muss be of class 'gn'")
  if (is.null(delta)) delta <- min(G$d)/2
  if (is.null(h)) h <- delta/2
  # line specific knot distances
  delta <- G$d*(delta > G$d) + delta*(delta <= G$d)
  delta <- pmin(G$d/floor(G$d/delta)*(G$d/delta - floor(G$d/delta) < 0.5) +
                  G$d/ceiling(G$d/delta)*(G$d/delta - floor(G$d/delta) >= 0.5),
                G$d/2)
  # ensure that h <= delta
  h <- min(h, min(delta))
  # line specific bin widths
  h <- G$d/floor(G$d/h)*(G$d/h - floor(G$d/h) < 0.5) +
    G$d/ceiling(G$d/h)*(G$d/h - floor(G$d/h) >= 0.5)
  # initializing...
  tau <- b <- z <- vector("list", G$M)
  N <- J <- rep(0, G$M)

  # do for every line segment
  for (m in 1:G$M) {
    # knot sequences tau
    tau[[m]] <- seq(0, G$d[m], delta[m])
    # bin boundaries b
    b[[m]] <- seq(0, G$d[m], h[m])
    # characterization of bins by midpoints z
    z[[m]] <- (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2
    # total count of bins in the geometric network
    N[m] <- length(z[[m]])
    # count of linear B-splines on line segment
    J[m] <- length(tau[[m]]) - 2
  }
  P <- list()
  P$splines <- list(delta = delta, tau = tau, J = J)
  P$bins <- list(h = h, b = b, z = z, N = N)
  P
}

#' Bin Point Pattern on a Geometric Network
#'
#' \code{bin_data} bins the data on the supplied point pattern according to
#' all possible combination of covariates.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param P as
#' @param vars as
#' @param intern as
#' @return The binned data.
#' @import dplyr
#' @importFrom stats setNames
#' @export

bin_data <- function(X, P, vars, intern){
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
  e <- NULL
  out <- list()
  if (is.element("x", vars)) {
    x <- rep(NA, sum(P$bins$N))
    ind <- 1
    test <- rep(NA, X$network$M)
    for (m in 1:X$network$M) {
      lins_m <- filter(X$network$lins, e == m)
      tp <- P$bins$z[[m]]/X$network$d[m]
      dx <- lins_m$v2_x - lins_m$v1_x
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        x[ind:(ind + length(ind_id) - 1)] <- lins_m$v1_x[i] + tp[ind_id]*dx[i]
        ind <- ind + length(ind_id)
      }
      test[m] <- ind - 1
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

