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
#' @importFrom stats setNames
#' @export

binData <- function(X, P, lins, smooths){
  # get covariates from every point on the network (if applicable)
  variables <- c(smooths, lins)
  if (all(variables %in% c(colnames(X$data)[-(1:5)], internal()))) {
    covariates <- as_tibble(X$data) %>% select(all_of(variables))
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
  dat <- tibble(id = rep(1:sum(P$bins$N), nrow(vars_comb)),
                count = NA, h = NA)
  dat <- bind_cols(dat, vars_comb %>% slice(rep(1:n(), each = sum(P$bins$N))))

  ind <- 1
  for (j in 1:nrow(vars_comb)) {
    ind_comb <- which(do.call(paste, covariates) == do.call(paste, vars_comb[j, ]))
    dat_sub <- X$data[ind_comb, ]
    for (m in 1:X$network$M) {
      # positions of data on curve e
      ind_e <- which(dat_sub$e == m)
      y_e <- sort(as.numeric(dat_sub[ind_e, ]$tp))*X$network$d[m]

      # bin data
      y_b <- rep(0, length(P$bins$z[[m]]))
      for (k in 1:length(P$bins$z[[m]])) {
        y_b[k] <- length(which(y_e < P$bins$b[[m]][k+1] & y_e > P$bins$b[[m]][k]))
      }

      # stack into one vector
      dat$count[ind:(ind + length(y_b) - 1)] <- y_b
      dat$h[ind:(ind + length(y_b) - 1)] <- P$bins$h[m]
      ind <- ind + length(y_b)
    }
    # add bin id for every row
    dat$id[((j-1)*sum(P$bins$N) + 1):(j*sum(P$bins$N))] <- 1:sum(P$bins$N)
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
    theta <- scoring(theta, design, rho)

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
#' @param offset as
#' @param delta The knot distance delta
#' @param h The bin width h
#' @param r The order of the penalty
#' @param density asd
#' @return an object of class gnppfit
#' @importFrom stats update as.formula model.matrix
#' @importFrom utils tail
#' @export

intensityPspline <- function(formula, X, offset = NULL,
                             delta = NULL, h = NULL, r = 1,
                             density = FALSE){
  # remove response from formula
  formula <- update(formula, NULL ~ .)

  # get linear and smooth terms
  variables <- all.vars(formula)
  if (length(variables) == 0){
    lins <- smooths <- NULL
  } else {
    variables2 <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
    lins <- variables[which(substring(variables2, 1, 2) != "s(")]
    smooths <- variables[which(substring(variables2, 1, 2) == "s(")]
  }
  # get representation of P-splines on the network
  P <- Pspline(X, delta, h, r)
  dat <- binData(X, P, lins, smooths)
  B <- K <- ind_smooths <- vector("list", length(smooths) + 1)
  K[[1]] <- P$K
  B[[1]] <- P$B

  # design matrix of the whole model
  Z <- P$B[dat$id, ]
  ind_smooths[[1]] <- 1:ncol(Z)
  names_theta <- paste0("G.", 1:ncol(Z))
  ind_lins <- NULL

  sm <- vector("list", length(smooths))
  if (length(sm) > 0) {
    for (i in 1:length(smooths)) {
      sm[[i]] <- smooth.construct(eval(parse(text = variables2[i])),
                                  data = dat, knots = NULL)
      D <- sm[[i]]$D[, -1]
      K[[i + 1]] <- t(D)%*%D
      Z_smooth <- sm[[i]]$X
      B[[i + 1]] <- sweep(Z_smooth, 2, colMeans(Z_smooth))[, -1]
      Z <- cbind(Z, sweep(Z_smooth, 2, colMeans(Z_smooth))[, -1])
      ind_smooths[[i + 1]] <- (ncol(Z) - ncol(Z_smooth) + 2):ncol(Z)
      names_theta <- c(names_theta, paste0(smooths[i], ".", 1:(ncol(Z_smooth) - 1)))
    }
  }

  # smooth effects
  # if (length(smooths) > 0){
  #   m <- 8
  #   for (a in 1:length(smooths)) {
  #     dat_a <- dat %>% pull(as.symbol(smooths[a]))
  #     z <- sort(unique(dat_a))
  #     delta <- (max(z) - min(z))/(m-1)
  #     k <- seq(min(z) - 3*delta, max(z) + 3*delta, delta)
  #     B.uncentered <- splineDesign(knots = k, x = z, ord = 4, outer.ok = TRUE)
  #
  #     #C <- rep(1, nrow(B.uncentered)) %*% B.uncentered
  #     #qrc <- qr(t(C))
  #     #Z1 <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]
  #     #B[[a + 1]] <- B.uncentered%*%Z1
  #
  #
  #     # remove first column for identification
  #     B[[a + 1]] <- sweep(B.uncentered, 2, colMeans(B.uncentered))[, -1]
  #     # compute penalty
  #     D <- diff(diag(ncol(B[[a + 1]]) + 1), differences = 2)[, -1]
  #     K[[a + 1]] <- t(D)%*%D
  #     # add to design matrix
  #     Z_test <- cbind(Z, B[[a + 1]][match(dat %>% pull(as.symbol(smooths[a])), z), ])
  #     ind_smooths[[a + 1]] <- (ncol(Z) - ncol(B[[a + 1]]) + 1):ncol(Z)
  #     names_theta <- c(names_theta, paste0(smooths[a], ".", 1:ncol(B[[a + 1]])))
  #   }
  # }

  if (length(lins) > 0){
    ff <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model.matrix.lin <- model.matrix(ff, data = dat)
    Z <- cbind(Z, model.matrix.lin[, -1])
    ind_lins <- (ncol(Z) - ncol(model.matrix.lin) + 2):ncol(Z)
    names_theta <- c(names_theta, colnames(model.matrix.lin)[-1])
  }

  dat$offset <- 1
  design <- list()
  design$K <- K
  design$data <- dat
  design$Z <- Z
  design$B <- B
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

  effects$smooth <- vector("list", length(smooths))
  names(effects$smooth) <- smooths
  if (length(smooths) > 0){
    for (i in 1:length(smooths)) {
      ind <- match(unique(design$data[[smooths[i]]]), design$data[[smooths[i]]])
      confidence_band <- getConfidenceBand(fit$theta, fit$V, design, i, ind, smooths)
      effects$smooth[[i]] <- tibble(x = design$data[[smooths[i]]][ind],
                                    y = as.vector(design$Z[, design$ind_smooths[[i + 1]]]%*%fit$theta[design$ind_smooths[[i + 1]]])[ind],
                                    lwr = confidence_band$lower,
                                    upr = confidence_band$upper)
    }
  }

  fit$effects <- effects
  fit$P <- P
  fit$data <- X$data
  fit$network <- X$network
  class(fit) <- "gnppfit"
  fit
}

#' Intensity Estimation on Geometric Networks with Penalized Splines
#'
#' \code{intensityPspline} estimates the intensity of a point pattern on a
#' geometric network.
#'
#' @param theta Point pattern on a geometric network (object of class \code{gnpp})
#' @param V Point pattern on a geometric network (object of class \code{gnpp})
#' @param design as
#' @param i as
#' @param ind as
#' @param smooths The knot distance delta
#' @param q The bin width h
#' @param R The order of the penalty
#' @return an object of class gnppfit
#' @importFrom stats quantile
#' @export

getConfidenceBand <- function(theta, V, design, i, ind, smooths, q = 0.05, R = 10000){
  gamma <- theta[design$ind_smooths[[i + 1]]]
  cov <- V[design$ind_smooths[[i + 1]], design$ind_smooths[[i + 1]]]
  x <- unique(design$data[[smooths[i]]])
  cov <- as.matrix(Matrix::Matrix(cov))

  set.seed(1)
  mu.sim <- matrix(0, R, length(x))
  for (j in 1:R) {
    gamma.sim <- mgcv::rmvn(1, gamma, cov)
    mu.sim[j, ] <- design$B[[i + 1]][ind, ]%*%gamma.sim
  }
  lower <- upper <- rep(0, ncol(mu.sim))
  for (j in 1:ncol(mu.sim)) {
    lower[j] <- stats::quantile(mu.sim[, j], probs = q/2)
    upper[j] <- stats::quantile(mu.sim[, j], probs = 1-q/2)
  }
  list(lower = lower, upper = upper)
}
