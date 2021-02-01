getDesign <- function(Gpp, smooths = NULL, lins = NULL, offset = NULL, m = 10, l = 2){

  Gpp <- addSpline(Gpp, delta = NULL, h = NULL, r = 1)

  K <- B <- ind_smooths <- vector("list", length(smooths) + 1)
  K[[1]] <- Gpp$K
  B[[1]] <- Gpp$B

  dat <- binData(Gpp, smooths = smooths, lins = lins)

  # baseline intensity of the network
  Z <- B[[1]][dat$id, ]
  ind_smooths[[1]] <- 1:ncol(Z)
  names_theta <- paste0("G.", 1:ncol(Z))

  # add network covariates if needed
  if (is.element("dist2V", c(smooths, lins))){
    data$dist2V <- dist2V(L)[data$id]
  }
  if (is.element("reldist2V", c(smooths, lins))){
    data$reldist2V <- dist2V(L, rel = TRUE)[data$id]
  }
  if (is.element("reldist2Vfrom", c(smooths, lins))){
    data$reldist2Vfrom <- tp(L)[data$id]
  }
  if (is.element("xcoord", c(smooths, lins))){
    data$xcoord <- get.x(L)[data$id]
  }
  if (is.element("x.km", c(smooths, lins))){
    data$x.km <- get.x.km(L)[data$id]
  }
  if (is.element("y.km", c(smooths, lins))){
    data$y.km <- get.y.km(L)[data$id]
  }
  if ("dist2Vdiscrete" %in% c(smooths, lins)){
    data$dist2Vdiscrete <- dist2Vdiscrete(L)[data$id]
  }
  if ("routetype" %in% lins){
    data$routetype <- L$routetype[data$id]
  }
  if ("direction" %in% lins){
    data$direction <- L$direction[data$id]
  }

  # smooth effects
  if (length(smooths) > 0){
    for (a in 1:length(smooths)) {
      data.a <- data %>% pull(as.symbol(smooths[a]))
      z <- sort(unique(data.a))
      delta <- (max(z) - min(z))/(m-1)
      k <- seq(min(z) - l*delta, max(z) + l*delta, delta)
      B.uncentered <- splineDesign(knots = k, x = z, ord = l+1, outer.ok = TRUE)
      # remove first column for identification
      B[[a + 1]] <- sweep(B.uncentered, 2, colMeans(B.uncentered))[, -1]
      # compute penalty
      D <- diff(diag(ncol(B[[a + 1]]) + 1), differences = 2)[, -1]
      K[[a + 1]] <- t(D)%*%D
      # add to design matrix
      Z <- cbind(Z, B[[a + 1]][match(data %>% pull(as.symbol(smooths[a])), z), ])
      ind.smooths[[a + 1]] <- (ncol(Z) - ncol(B[[a + 1]]) + 1):ncol(Z)
      names.theta <- c(names.theta, paste0(smooths[a], ".", 1:ncol(B[[a + 1]])))
    }
  }

  # linear effects
  ind_lins <- NULL
  if (length(lins) > 0){
    ff <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model.matrix.lin <- model.matrix(ff, data = data)
    Z <- cbind(Z, model.matrix.lin[, -1])
    ind.lins <- (ncol(Z) - ncol(model.matrix.lin) + 2):ncol(Z)
    names.theta <- c(names.theta, colnames(model.matrix.lin)[-1])
  }

  # offset
  if (!is.null(offset)){
    dat$offset <- get.offset(offset, L)[data$id]
  } else {
    dat$offset <- 1
  }

  list(data = dat, Z = Z, B = B, K = K, ind_smooths = ind_smooths, ind_lins = ind_lins, names_theta = names_theta)
}

binData <- function(Gpp, smooths = NULL, lins = NULL){

  # get covariates from every point on the network (if applicable)
  name <- intersect(names(Gpp$data), c(smooths, lins))
  covariates <- as.data.frame(Gpp$data) %>% dplyr::select(dplyr::all_of(name))
  if (ncol(covariates) == 0){
    covariates <- covariates %>% dplyr::mutate(q = 1)
  }

  # get all combinations of covariates and calculate the number of rows of the data matrix
  covariates.comb <- covariates %>%
    dplyr::distinct() %>%
    expand.grid() %>%
    dplyr::distinct() %>%
    dplyr::as_tibble()

  for (a in 1:length(covariates.comb)) {
    covariates.comb <- dplyr::arrange(covariates.comb, !!dplyr::sym(names(covariates.comb[a])))
  }

  N <- sum(Gpp$bins$N)*nrow(covariates.comb)

  # initializing
  dat <- stats::setNames(data.frame(matrix(nrow = N, ncol = length(name) + 3)), c("id", "count", "h", name)) %>% dplyr::as_tibble() %>%
    dplyr::mutate(id = as.integer(id), count = as.double(count), h = as.double(h))

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
    data.sub <- Gpp$data[ind.cov, ]
    for (m in 1:Gpp$M) {
      # positions of data on line m
      ind.m <- which(data.sub$e == m)
      y.m <- sort(as.numeric(data.sub[ind.m, ]$tp))*Gpp$d[m]

      # bin data
      y.b <- rep(0, length(Gpp$bins$z[[m]]))
      for (k in 1:length(Gpp$bins$z[[m]])) {
        y.b[k] <- length(which(y.m < Gpp$bins$b[[m]][k+1] & y.m > Gpp$bins$b[[m]][k]))
      }

      # stack into one vector
      dat$count[ind:(ind + length(y.b) - 1)] <- y.b
      dat$h[ind:(ind + length(y.b) - 1)] <- Gpp$bins$h[m]
      ind <- ind + length(y.b)
    }
    # add bin id for every row
    dat$id[((j-1)*sum(Gpp$bins$N) + 1):(j*sum(Gpp$bins$N))] <- 1:sum(Gpp$bins$N)

    # add covariates
    if (ncol(dat) > 3){
      dat[((j-1)*sum(Gpp$bins$N) + 1):(j*sum(Gpp$bins$N)), 4:ncol(dat)] <- covariates.comb[j, ]
    }
  }
  dat
}

# fit.lpp = function(L.lpp, smooths = NULL, lins = NULL, offset = NULL,
#                    rho = 10, rho.max = 1e5, eps.rho = 0.01, maxit.rho = 100){
fitData <- function(G, smooths = NULL, lins = NULL, offset = NULL, rho = 10,
                    rho_max = 1e5, eps_rho = 0.01, maxit_rho = 100){
  design <- getDesign(G)

  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(rho, length(design$K))
  Delta_rho <- Inf
  it_rho <- 0
  theta <- rep(0, ncol(design$Z))
  while(Delta_rho > eps_rho){
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
    print(rho_new)
    Delta_rho <- sqrt(sum((rho_new - rho)^2))/sqrt(sum((rho)^2))
    rho <- rho_new
  }

  # effects in one table
  effects <- list(linear = NULL, smooth = NULL)

  #linear effects
  if (length(lins) > 0) {
    effects$linear <- tibble(name = tail(design$names.theta, length(design$ind.lins)), estimate = NA, se = NA)
    for (i in 1:length(design$ind.lins)) {
      effects$linear$estimate[i] <- round(theta[design$ind.lins[i]], 3)
      effects$linear$se[i] <- round(sqrt(V[design$ind.lins[i], design$ind.lins[i]]), 3)
      effects$linear$rr[i] <- round(exp(effects$linear$estimate[i]), 2)
      effects$linear$rr.lower[i] <- round(exp(effects$linear$estimate[i] - 1.96*effects$linear$se[i]), 2)
      effects$linear$rr.upper[i] <- round(exp(effects$linear$estimate[i] + 1.96*effects$linear$se[i]), 2)
    }
  }

  # smooth effects
  effects$smooth <- vector("list", length(smooths))
  names(effects$smooth) <- smooths
  if (length(smooths) > 0){
    for (i in 1:length(smooths)) {
      ind <- match(unique(design$data[[smooths[i]]]), design$data[[smooths[i]]])
      confidence.band <- get.confidence.band(theta, V, design, i, ind, smooths)
      effects$smooth[[i]] <- tibble(x = design$data[[smooths[i]]][ind],
                                    y = as.vector(design$Z[, design$ind.smooths[[i + 1]]]%*%theta[design$ind.smooths[[i + 1]]])[ind],
                                    lwr = confidence.band$lower,
                                    upr = confidence.band$upper)
    }
  }

  list(theta = theta, V = V, ind_smooths = design$ind_smooths, names_theta = design$names_theta, effects = effects)
}

intensity <- function(G, density = FALSE){
  fit <- fitData(G)
  fit
}
