internal <- function(){
  c("dist2V")
}

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
  # get covariates from every point on the network (if applicable)
  vars <- c(smooths, lins)
  if (all(vars %in% c(colnames(X$data)[-(1:5)], internal()))) {
    covariates <- as_tibble(X$data) %>% select(all_of(vars))
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
  dat <- tibble(seg = rep(1:sum(P$bins$N), nrow(vars_comb)),
                count = NA, h = NA)
  dat <- bind_cols(dat, vars_comb %>% slice(rep(1:n(), each = sum(P$bins$N))))

  ind <- 1
  for (j in 1:nrow(vars_comb)) {
    ind_comb <- which(do.call(paste, list(vars)) == do.call(paste, vars_comb[j, ]))
    data.sub <- X$data[ind_comb, ]
    for (m in 1:X$network$M) {
      # positions of data on curve e
      ind_e <- which(data.sub$e == m)
      y.m <- sort(as.numeric(data.sub[ind_e, ]$tp))*X$network$d[m]

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
    dat$seg[((j-1)*sum(P$bins$N) + 1):(j*sum(P$bins$N))] <- 1:sum(P$bins$N)
  }
  dat
}

#' Fit a Penalized Spline Poisson Model on a Geometric Network
#'
#' \code{fitModel} performs the actual model fit.
#'
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

fitModel <- function(design, lins = NULL, smooths = NULL, offset = NULL, rho = 10,
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
    V <- Matrix::solve(fisher(theta, design, rho))

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
#' @param formula Point pattern on a geometric network (object of class \code{gnpp})
#' @param X Point pattern on a geometric network (object of class \code{gnpp})
#' @param lins as
#' @param smooths as
#' @param offset as
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @param density asd
#' @return an object of class gnppfit
#' @importFrom stats update
#' @export

intensityPspline <- function(formula, X, offset = NULL,
                             delta = NULL, h = NULL, r = 1,
                             density = FALSE){
  # remove response from formula
  formula <- update(formula, NULL ~ .)

  # get linear and smooth variables
  vars <- all.vars(formula)
  vars2 <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
  lins <- vars[which(substring(vars2, 1, 2) != "s(")]
  smooths <- vars[which(substring(vars2, 1, 2) == "s(")]

  P <- Pspline(X, delta, h, r)
  K <- B <- ind_smooths <- vector("list", length(smooths) + 1)
  design <- list()
  K[[1]] <- P$K
  B[[1]] <- P$B
  dat <- binData(X, P, lins, smooths)

  Z <- B[[1]][dat$seg, ]
  ind_smooths[[1]] <- 1:ncol(Z)
  names_theta <- paste0("G.", 1:ncol(Z))
  ind_lins <- NULL

  if (length(lins) > 0){
    ff <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model.matrix.lin <- model.matrix(ff, data = dat)
    Z <- cbind(Z, model.matrix.lin[, -1])
    ind_lins <- (ncol(Z) - ncol(model.matrix.lin) + 2):ncol(Z)
    names_theta <- c(names_theta, colnames(model.matrix.lin)[-1])
  }

  dat$offset <- 1

  design$K <- K
  design$B <- B
  design$data <- dat
  design$Z <- Z
  design$ind_smooths <- ind_smooths
  design$names_theta <- names_theta

  fit <- fitModel(design, lins, smooths, offset)

  # effects in one table
  effects <- list(linear = NULL, smooth = NULL)

  #linear effects
  if (length(lins) > 0) {
    effects$linear <- tibble(name = tail(design$names_theta, length(ind_lins)),
                             estimate = NA, se = NA, rr = NA, rr.lower = NA, rr.upper = NA)
    for (i in 1:length(ind_lins)) {
      effects$linear$estimate[i] <- round(fit$theta[ind_lins[i]], 3)
      effects$linear$se[i] <- round(sqrt(fit$V[ind_lins[i], ind_lins[i]]), 3)
      effects$linear$rr[i] <- round(exp(effects$linear$estimate[i]), 2)
      effects$linear$rr.lower[i] <- round(exp(effects$linear$estimate[i] - 1.96*effects$linear$se[i]), 2)
      effects$linear$rr.upper[i] <- round(exp(effects$linear$estimate[i] + 1.96*effects$linear$se[i]), 2)
    }
  }
  fit$effects <- effects
  fit$P <- P
  fit$data <- X$data
  fit$network <- X$network
  class(fit) <- "gnppfit"
  fit
}
