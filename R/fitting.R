#' Bin Point Pattern on a Geometric Network
#'
#' \code{binData} performs the actual model fit.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param P as
#' @param lins as
#' @param smooths as
#' @return The binned data.
#' @import dplyr
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @export

binData <- function(X, P, lins, smooths){
  h <- NULL

  # get covariates from every point on the network (if applicable)
  name <- intersect(names(X$data), c(smooths, lins))
  covariates <- as.data.frame(X$data) %>% select(all_of(name))
  if (ncol(covariates) == 0){
    covariates <- covariates %>%mutate(q = 1)
  }

  # get all combinations of covariates and calculate the number of rows of the data matrix
  covariates.comb <- covariates %>%
    distinct() %>%
    expand.grid() %>%
    distinct() %>%
    as_tibble()

  for (a in 1:length(covariates.comb)) {
    covariates.comb <- arrange(covariates.comb, !!sym(names(covariates.comb[a])))
  }

  N <- sum(P$bins$N)*nrow(covariates.comb)

  # initializing
  dat <- setNames(data.frame(matrix(nrow = N, ncol = length(name) + 3)), c("id", "count", "h", name)) %>% as_tibble() %>%
    mutate(id = as.integer(id), count = as.double(count), h = as.double(h))

  # set factor variable if applicable
  if (length(name) > 0){
    for (a in 1:length(name)) {
      if (is.factor(covariates.comb %>% pull(sym(name[a])))){
        dat <- mutate(dat, !!name[a] := factor(NA, levels = levels(covariates.comb %>% pull(sym(name[a])))))
      }
      if (is.double(covariates.comb %>% pull(sym(name[a])))){
        dat <- mutate(dat, !!name[a] := as.double(NA))
      }
      if (is.integer(covariates.comb %>% pull(sym(name[a])))){
        dat <- mutate(dat, !!name[a] := as.integer(NA))
      }
    }
  }

  ind <- 1
  for (j in 1:nrow(covariates.comb)) {
    ind.cov <- which(do.call(paste, covariates) == do.call(paste, covariates.comb[j, ]))
    data.sub <- X$data[ind.cov, ]
    for (m in 1:X$M) {
      # positions of data on line m
      ind.m <- which(data.sub$e == m)
      y.m <- sort(as.numeric(data.sub[ind.m, ]$tp))*X$d[m]

      # bin data
      y.b <- rep(0, length(P$bins$z[[m]]))
      for (k in 1:length(P$bins$z[[m]])) {
        y.b[k] <- length(which(y.m < P$bins$b[[m]][k+1] & y.m > P$bins$b[[m]][k]))
      }

      # stack into one vector
      dat$count[ind:(ind + length(y.b) - 1)] <- y.b
      dat$h[ind:(ind + length(y.b) - 1)] <- P$bins$h[m]
      ind <- ind + length(y.b)
    }
    # add bin id for every row
    dat$id[((j-1)*sum(P$bins$N) + 1):(j*sum(P$bins$N))] <- 1:sum(P$bins$N)

    # add covariates
    if (ncol(dat) > 3){
      dat[((j-1)*sum(P$bins$N) + 1):(j*sum(P$bins$N)), 4:ncol(dat)] <- covariates.comb[j, ]
    }
  }
  dat
}

#' Fit a Penalized Spline Poisson Model on a Geometric Network
#'
#' \code{fitModel} performs the actual model fit.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param design as
#' @param lins as
#' @param smooths as
#' @param offset as
#' @param rho as
#' @param rho_max as
#' @param eps_rho as
#' @param maxit_rho as
#' @return Model fit.
#' @import dplyr
#' @importFrom stats optim
#' @export

fitModel <- function(X, design, lins = NULL, smooths = NULL, offset = NULL, rho = 10,
                    rho_max = 1e5, eps_rho = 0.01, maxit_rho = 100){

  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(rho, length(design$K))
  Delta_rho <- Inf
  it_rho <- 0
  theta <- rep(0, ncol(design$Z))
  while(Delta_rho > eps_rho){
    print(rho)
    it_rho <- it_rho + 1
    fit <- optim(theta, fn = logL, gr = score, design = design, rho = rho,
                 control = list(fnscale = -1, maxit = 1000, factr = 1e4),
                 method = "L-BFGS-B")
    theta <- fit$par
    V <- solve(fisher(theta, design, rho))

    # update rho
    rho_new <- rep(NA, length(design$K))
    for (a in 1:length(design$K)) {
      rho_new[a] <- as.vector(rho[a]*(Matrix::rankMatrix(design$K[[a]], method = "qr.R")[1]/rho[a] -
                                        sum(Matrix::diag(V[design$ind_smooths[[a]], design$ind_smooths[[a]]]%*%design$K[[a]])))/
                                (Matrix::t(theta[design$ind_smooths[[a]]])%*%design$K[[a]]%*%theta[design$ind_smooths[[a]]]))
    }

    if (any(rho_new > rho_max)) break
    if (any(rho_new < 0)){
      warning("rho = 0 has occurred")
    }
    if (it_rho > maxit_rho){
      warning("Stopped estimation of rho because maximum number of iterations has been reached!")
      break
    }
    Delta_rho <- sqrt(sum((rho_new - rho)^2))/sqrt(sum((rho)^2))
    rho <- rho_new
  }

  # effects in one table
  effects <- list(linear = NULL, smooth = NULL)

  #linear effects
  if (length(lins) > 0) {
    effects$linear <- tibble(name = utils::tail(design$names.theta, length(design$ind.lins)), estimate = NA, se = NA)
  }

  # smooth effects
  effects$smooth <- vector("list", length(smooths))
  names(effects$smooth) <- smooths

  list(theta = theta, V = V, ind_smooths = design$ind_smooths, names_theta = design$names_theta, effects = effects)
}

#' Intensity Estimation on Geometric Networks with Penalized Splines
#'
#' \code{intensityPspline} estimates the intensity of a point pattern on a
#' geometric network.
#'
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param lins as
#' @param smooths as
#' @param offset as
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @param density asd
#' @return an object of class gnppfit
#' @export

intensityPspline <- function(X, lins = NULL, smooths = NULL, offset = NULL,
                             delta = NULL, h = NULL, r = 1, density = FALSE){

  P <- Pspline(X, delta, h, r)
  K <- B <- ind_smooths <- vector("list", length(smooths) + 1)
  design <- list()
  K[[1]] <- P$K
  B[[1]] <- P$B
  dat <- binData(X, P, lins, smooths)

  Z <- B[[1]][dat$id, ]
  ind_smooths[[1]] <- 1:ncol(Z)
  names_theta <- paste0("G.", 1:ncol(Z))
  ind_lins <- NULL

  dat$offset <- 1

  design$K <- K
  design$B <- B
  design$data <- dat
  design$Z <- Z
  design$ind_smooths <- ind_smooths
  design$names_theta <- names_theta

  fit <- fitModel(X, design, lins, smooths, offset)
  X$model <- list(fit = fit, P = P)
  class(X) <- "gnppfit"
  X
}
