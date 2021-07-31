#' Computes a Global Knot Distance from the Input
#'
#' @param G An object of class \code{gn}) or a point pattern on a geometric
#' network (object of class \code{gnpp}).
#' @param delta A numeric vector of length one which already defines
#' the global knot distance. Alternatively, delta can be supplied in terms of a
#' quantile of the curve lengths of the network, i.e. a number in the unit
#' interval. In the latter case, delta must be supplied as a character vector of
#' length one, see the examples. By default, delta is chosen to be half of the
#' minimal curve length.
#' @param h A numeric vector of length one which already defines the
#' global knot distance h. Alternatively, h can be supplied in terms of a
#' fraction of delta, i.e. a number in the interval (0,1]. In the latter case,
#' h must be supplied as character vector of length one, see the examples.
#' By default, h is chosen to be half of the global knot distance.
#' @return The global knot distance delta.
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export
#' @examples
#' G <- as_gn(montgomery)
#' # use default arguments
#' setup <- delta_h_global(G)
#' setup
#' # set numeric value for delta and fraction for h
#' setup <- delta_h_global(G, delta = 0.1, h = "0.25")
#' setup
#' # set quantile for delta
#' setup <- delta_h_global(G, delta = "0.05")
#' setup

delta_h_global <- function(G, delta = NULL, h = NULL) {
  #global knot distance
  if (is.null(delta)) {
    delta <- min(G$d)/2
  } else if (is.character(delta) & !is.na(as.numeric(delta))) {
    if (is.character(delta) & !is.na(as.numeric(delta))) {
      quant <- as.numeric(delta)
      if (quant >= 0 & quant <= 1) {
        delta <- as.numeric(quantile(G$d, quant)/2)
        if (delta > quantile(G$d, 0.5)/2) {
          warning(paste0("Many network segments have no curve specific ",
                         "B-splines.\nChoose a lower value for the quantile."))
        }
      } else {
        stop("The quantile supplied for delta must be in the interval [0,1]")
      }
    } else {
      stop("The character must be convertible to class numeric.")
    }
  } else if (is.numeric(delta)) {
    if (delta <= 0) {
      stop("delta must be positive!")
    }
    if (delta > quantile(G$d, 0.5)/2) {
      warning(paste0("Many network segments have no curve specific ",
                     "B-splines.\nChoose a lower value for delta."))
    }
  } else {
    stop("Supply delta in correct format, see the documentation.")
  }
  if (is.null(h)) {
    h <- delta/2
  } else if (is.numeric(h)) {
    if (h <= 0) stop("h must be positive!")
    if (h > delta) stop("h must be at most as large as delta!")
  } else if (is.character(h) & !is.na(as.numeric(h))) {
    mult <- as.numeric(h)
    h <- delta*mult
  } else {
    stop("Supply h in correct format, see the documentation.")
  }
  list(delta = delta, h = h)
}

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

network_knots <- function(G, delta){
  if (inherits(G, "gnpp")) G <- as_gn(G)
  if (!inherits(G, "gn")) stop("G muss be of class 'gn'")
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
  # line specific bin widths
  h <- pmin(G$d/floor(G$d/h)*(G$d/h - floor(G$d/h) < 0.5) +
    G$d/ceiling(G$d/h)*(G$d/h - floor(G$d/h) >= 0.5), G$d/2)
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
      y_e <- sort(as.numeric(data_sub[ind_e, ]$tp_e))*X$network$d[m]

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
      dx <- lins_m$a2_x - lins_m$a1_x
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        x[ind:(ind + length(ind_id) - 1)] <- lins_m$a1_x[i] + tp[ind_id]*dx[i]
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
      dy <- lins_m$a2_y - lins_m$a1_y
      for (i in 1:nrow(lins_m)) {
        ind_id <- which(tp >= lins_m$frac1[i] & tp < sum(lins_m$frac2[1:i]))
        y[ind:(ind + length(ind_id) - 1)] <- lins_m$a1_y[i] + tp[ind_id]*dy[i]
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
    if (!is.null(scale$dist2V)) out$dist2V <- d*scale$dist2V
    else out$dist2V <- d
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

#' Fitted Intensity on a Geometric Network
#'
#' @param z The shortest path distance from the beginning of the network
#' segment.
#' @param m The network segment index.
#' @param fit1 A fitted geometric network.
#' @param fit2 A second fitted geometric network. If specified, the function
#' returns the squared difference of the intensity fits at the specified point
#' of the network.
#' @param scale A numeric vector of length two which determines the scaling
#' of the two intensity functions.
#'
#' @return A numeric vector of length one, indicating the intensity (or the
#' squared difference of two intensities) at the specified point.
#' @export

network_intensity <- function(z, m, fit1, fit2 = NULL, scale = NULL){

  ind_v1 <- which(fit1$network$incidence[, m] == -1)
  ind_v2 <- which(fit1$network$incidence[, m] == 1)

  B1 <- matrix(0, length(z), fit1$knots$J[m] + 2)
  B1[, 1] <- (1-z/fit1$knots$delta[m])*(1-z/fit1$knots$delta[m] > 0)
  B1[, 2:(ncol(B1)-1)] <- splineDesign(knots = fit1$knots$tau[[m]], x = z, ord = 2, outer.ok = TRUE)
  B1[, ncol(B1)] <- (1-(fit1$network$d[m]-z)/fit1$knots$delta[m])*(1-(fit1$network$d[m]-z)/fit1$knots$delta[m] > 0)

  gamma1 <- fit1$coefficients[fit1$ind[[1]]]


  gamma1_m <- c(gamma1[sum(fit1$knots$J) + ind_v1],
                gamma1[((cumsum(fit1$knots$J)-fit1$knots$J)[m]+1):cumsum(fit1$knots$J)[m]],
                gamma1[sum(fit1$knots$J) + ind_v2])

  if (!is.null(fit2)) {
    B2 <- matrix(0, length(z), fit2$knots$J[m] + 2)
    B2[, 1] <- (1-z/fit2$knots$delta[m])*(1-z/fit2$knots$delta[m] > 0)
    B2[, 2:(ncol(B2)-1)] <- splineDesign(knots = fit2$knots$tau[[m]], x = z, ord = 2, outer.ok = TRUE)
    B2[, ncol(B2)] <- (1-(fit2$network$d[m]-z)/fit2$knots$delta[m])*(1-(fit2$network$d[m]-z)/fit2$knots$delta[m] > 0)

    gamma2 <- fit2$coefficients[fit2$ind[[1]]]

    gamma2_m <- c(gamma2[sum(fit2$knots$J) + ind_v1],
                  gamma2[((cumsum(fit2$knots$J)-fit2$knots$J)[m]+1):cumsum(fit2$knots$J)[m]],
                  gamma2[sum(fit2$knots$J) + ind_v2])
  }

  if (min(z) >= 0 & max(z) <= fit1$network$d[m]) {
    if (!is.null(fit2)) {
      if (is.null(scale)) scale = c(1, 1)
      return(as.vector((exp(B1%*%gamma1_m)/scale[1] - exp(B2%*%gamma2_m)/scale[2])^2))
    } else {
      return(as.vector(exp(B1%*%gamma1_m)))
    }
  } else {
    stop("Wrong domain of z on this edge!")
  }
}

#' Integral of a fitted intensity
#'
#' @param fit A fitted point process on a geometric network.
#'
#' @return The integral.
#' @export

network_integral <- function(fit) {
  int <- rep(NA, fit$network$M)
  for (m in 1:fit$network$M) {
    int[m] <- integrate(network_intensity, lower = 0, upper = fit$network$d[m],
                        m = m, fit1 = fit)$value
  }
  sum(int)
}

#' Computation of the Integrated Squared Error
#'
#' Computes the integrated squared error between a true between a true
#' intensity \code{fit1} and an estimate \code{fit2}.
#'
#' @param fit1 The true intensity.
#' @param fit2 The estimated intensity.
#'
#' @return The (normalized) integrated squared error.
#' @export

network_ISE <- function(fit1, fit2) {

  scale <- c(nrow(fit1$data), network_integral(fit2))
  int <- rep(NA, fit1$network$M)
  for (m in 1:fit1$network$M) {
    int[m] <- integrate(network_intensity, lower = 0, upper = fit1$network$d[m],
                        m = m, fit1 = fit1, fit2 = fit2, scale = scale)$value
  }
  sum(int)
}

#' Find Location of a Point on a Geometric Network
#'
#' @param G A geometric network.
#' @param m The segment index.
#' @param z The shortest path distance from the beginning of the network
#' segment.
#'
#' @return A list with names x and y which contains the coordinated of the
#' points.
#' @export

network_location <- function(G, m, z){
  e <- NULL
  lins_m <- filter(G$lins, e == m)
  frac <- c(lins_m$frac1, 1)
  tp_e <- z/G$d[m]
  lin <- findInterval(tp_e, frac, rightmost.closed = TRUE, all.inside = TRUE)
  l <- lins_m$l[lin]
  dx <- lins_m$a2_x[lin] - lins_m$a1_x[lin]
  dy <- lins_m$a2_y[lin] - lins_m$a1_y[lin]
  tp_l <- (z - lins_m$frac1[lin]*G$d[m])/(lins_m$length[lin])
  x <- lins_m$a1_x[lin] + tp_l*dx
  y <- lins_m$a1_y[lin] + tp_l*dy
  list(l = l, tp_l = tp_l, x = x, y = y)
}



