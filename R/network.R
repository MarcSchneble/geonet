#' Defining knots on a Geometric Network
#'
#' \code{network_knots} defines knots on a geometric network (object of class
#' \code{gn}) which can be used to construct linear B-splines on it.
#'
#' @param G An object of class \code{gn}) or a point pattern on a geometric
#' network (object of class \code{gnpp}).
#' @param delta The global knot distance delta.
#' @return A list which contains the knot sequence of every curve of the
#' geometric network.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

network_knots <- function(G, delta = NULL){
  if (inherits(G, "gnpp")) G <- as_gn(G)
  if (!inherits(G, "gn")) stop("G muss be of class 'gn'")
  if (is.null(delta)) delta <- min(G$d)/2
  # line specific knot distances
  delta <- G$d*(delta > G$d) + delta*(delta <= G$d)
  delta <- pmin(G$d/floor(G$d/delta)*(G$d/delta - floor(G$d/delta) < 0.5) +
                  G$d/ceiling(G$d/delta)*(G$d/delta - floor(G$d/delta) >= 0.5),
                G$d/2)
  # initializing...
  tau <- vector("list", G$M)
  J <- rep(0, G$M)

  # do for every line segment
  for (m in 1:G$M) {
    # knot sequences tau
    tau[[m]] <- seq(0, G$d[m], delta[m])
    # count of linear B-splines on line segment
    J[m] <- length(tau[[m]]) - 2
  }
  knots <- list(delta = delta, tau = tau, J = J)
  knots
}

#' Defining bins on a Geometric Network
#'
#' \code{network_bins} subdivides each curve segment into several bins.
#'
#' @param G An object of class \code{gn}) or a point pattern on a geometric
#' network (object of class \code{gnpp}).
#' @param h The global bin width h.
#' @return A list which contains the bins of every curve of the geometric
#' network.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

network_bins <- function(G, h = NULL){
  if (inherits(G, "gnpp")) G <- as_gn(G)
  if (!inherits(G, "gn")) stop("G muss be of class 'gn'")
  if (is.null(h)) h <- min(G$d)/4
  # line specific bin widths
  h <- G$d/floor(G$d/h)*(G$d/h - floor(G$d/h) < 0.5) +
    G$d/ceiling(G$d/h)*(G$d/h - floor(G$d/h) >= 0.5)
  # initializing...
  b <- z <- vector("list", G$M)
  N <- rep(0, G$M)

  # do for every line segment
  for (m in 1:G$M) {
    # bin boundaries b
    b[[m]] <- seq(0, G$d[m], h[m])
    # characterization of bins by midpoints z
    z[[m]] <- (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2
    # total count of bins in the geometric network
    N[m] <- length(z[[m]])
  }
  bins <- list(h = h, b = b, z = z, N = N)
  bins
}

#' Bin Point Pattern on a Geometric Network
#'
#' \code{bin_data} bins the data on the supplied point pattern according to
#' all possible combination of covariates.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param bins A list containing the bins of the geometric network.
#' @param vars A character vector containing the name of all covariates in the
#' model.
#' @param vars_internal A character vector containing the name of all interval
#' covariates in the model.
#' @param scale A named list which specifies the rescaling of network related
#' covariates. Currently, only x- and y-coordinates can be scaled.
#' @return The binned data.
#' @import dplyr
#' @importFrom stats setNames
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

bin_data <- function(X, bins = NULL, vars = NULL, vars_internal = NULL,
                     scale = NULL){
  # get covariates from every point on the network (if applicable)
  if (is.null(bins)) bins <- network_bins(X$network)
  vars_external <- setdiff(vars, vars_internal)
  if (length(vars_external) > 0) {
    covariates <- as_tibble(X$data) %>% select(all_of(vars_external))
  } else {
    covariates <- tibble(a = rep(1, nrow(X$data)))
  }

  # get all combinations of covariates and calculate the number of rows of the data matrix
  vars_comb <- covariates %>% distinct() %>% expand.grid() %>%
    distinct() %>% as_tibble()

  # sort data frame
  for (a in 1:length(vars_comb)) {
    vars_comb <- arrange(vars_comb, !!sym(names(vars_comb[a])))
  }

  # initializing
  data <- tibble(id = 1:sum(bins$N), count = NA, h = NA)
  if (length(vars_internal) > 0) {
    data <- bind_cols(data, as_tibble(internal(vars_internal, X, bins, scale)))
  }
  data <- data %>% slice(rep(1:n(), nrow(vars_comb))) %>%
    bind_cols(vars_comb %>% slice(rep(1:n(), each = sum(bins$N))))

  ind <- 1
  for (j in 1:nrow(vars_comb)) {
    ind_comb <- which(do.call(paste, covariates) == do.call(paste, vars_comb[j, ]))
    data_sub <- X$data[ind_comb, ]
    for (m in 1:X$network$M) {
      # positions of data on curve e
      ind_e <- which(data_sub$e == m)
      y_e <- sort(as.numeric(data_sub[ind_e, ]$tp))*X$network$d[m]

      # bin data
      y_b <- rep(0, length(bins$z[[m]]))
      for (k in 1:length(bins$z[[m]])) {
        y_b[k] <- length(which(y_e < bins$b[[m]][k+1] & y_e > bins$b[[m]][k]))
      }

      # stack into one vector
      data$count[ind:(ind + length(y_b) - 1)] <- y_b
      data$h[ind:(ind + length(y_b) - 1)] <- bins$h[m]
      ind <- ind + length(y_b)
    }
    # add bin id for every row
    data$id[((j-1)*sum(bins$N) + 1):(j*sum(bins$N))] <- 1:sum(bins$N)
  }
  data$offset <- 1
  data
}


#' Internal Covariates
#'
#' \code{internal} computes the values of internal covariates at the midpoints
#' of the bins of the network. Internal covariates can either be supplied via
#' the point pattern or they are a function of the network. Currently, x- and y-
#' coordinates are supported for the latter.
#'
#' @param vars The name of the covariates which should go into the model
#' as linear internal covariates.
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param bins A list containing the bins of the geometric network.
#' @param scale A named list which specifies the rescaling of network related
#' covariates. Currently, only x- and y-coordinates can be scaled.
#' @return A data frame with the number of rows equal to the number of bins of the
#' geometric network (\code{sum(bins$N)}) and the number of columns equal to
#' the length of \code{vars}.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

internal <- function(vars, X, bins, scale){
  e <- NULL
  out <- list()
  if (is.element("x", vars)) {
    x <- rep(NA, sum(bins$N))
    ind <- 1
    test <- rep(NA, X$network$M)
    for (m in 1:X$network$M) {
      lins_m <- filter(X$network$lins, e == m)
      tp <- bins$z[[m]]/X$network$d[m]
      dx <- lins_m$v2_x - lins_m$v1_x
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        x[ind:(ind + length(ind_id) - 1)] <- lins_m$v1_x[i] + tp[ind_id]*dx[i]
        ind <- ind + length(ind_id)
      }
      test[m] <- ind - 1
    }
    if (!is.null(scale$x)) out$x <- x*scale$x
    else out$x <- x
  }
  if (is.element("y", vars)) {
    y <- rep(NA, sum(bins$N))
    ind <- 1
    for (m in 1:X$network$M) {
      lins_m <- filter(X$network$lins, e == m)
      tp <- bins$z[[m]]/X$network$d[m]
      dy <- lins_m$v2_y - lins_m$v1_y
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        y[ind:(ind + length(ind_id) - 1)] <- lins_m$v1_y[i] + tp[ind_id]*dy[i]
        ind <- ind + length(ind_id)
      }
    }
    if (!is.null(scale$y)) out$y <- y*scale$y
    else out$y <- y
  }
  if (is.element("dist2V", vars)) {
    d <- rep(NA, sum(bins$N))
    ind <- 1
    for (m in 1:X$network$M) {
      d[ind:(ind + bins$N[m] - 1)] <- pmin(bins$z[[m]], X$network$d[m] - bins$z[[m]])
      ind <- ind + bins$N[m]
    }
    out$dist2V <- d
  }
  X_internal <- colnames(X$network$lins[-(1:11)])
  if (length(X_internal) > 0) {
    for (k in 1:length(X_internal)) {
      if (is.element(X_internal[k], vars)) {
        if (class(X$network$lins[[X_internal[k]]]) == "factor") {
          val <- factor(rep(NA, sum(bins$N)),
                        levels = levels(X$network$lins[[X_internal[k]]]))
        } else {
          val <- rep(NA, sum(bins$N))
        }
        ind <- 1
        for (m in 1:X$network$M) {
          lins_m <- filter(X$network$lins, e == m)
          if (length(unique(lins_m[[paste(X_internal[k])]])) > 1) {
            #stop("Internal covariates must be unique on each curve.")
          }
          val[ind:(ind + bins$N[m] - 1)] <- lins_m[[paste(X_internal[k])]][1]
          ind <- ind + bins$N[m]
        }
        out[[paste(X_internal[k])]] <- val
      }
    }
  }
  out
}

