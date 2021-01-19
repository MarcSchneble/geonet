getDesign <- function(Gpp, smooths = NULL, lins = NULL, offset = NULL, m = 10, l = 2){

  Gpp <- addSpline(Gpp, delta = NULL, h = NULL, r = 1)

  K <- B <- ind_smooths <- vector("list", length(smooths) + 1)
  K[[1]] <- Gpp$K
  B[[1]] <- Gpp$B

  dat <- binData(Gpp, smooths = smooths, lins = lins)

  # baseline intensity of the network
  Z <- B[[1]][data$id, ]
  ind.smooths[[1]] <- 1:ncol(Z)
  names.theta <- paste0("L.", 1:ncol(Z))

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
  ind.lins <- NULL
  if (length(lins) > 0){
    ff <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model.matrix.lin <- model.matrix(ff, data = data)
    Z <- cbind(Z, model.matrix.lin[, -1])
    ind.lins <- (ncol(Z) - ncol(model.matrix.lin) + 2):ncol(Z)
    names.theta <- c(names.theta, colnames(model.matrix.lin)[-1])
  }

  # offset
  if (!is.null(offset)){
    data$offset <- get.offset(offset, L)[data$id]
  } else {
    data$offset <- 1
  }

  return(list(data = data, Z = Z, B = B, K = K, ind.smooths = ind.smooths, ind.lins = ind.lins, names.theta = names.theta))
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
  data <- stats::setNames(data.frame(matrix(nrow = N, ncol = length(name) + 3)), c("id", "count", "h", name)) %>% dplyr::as_tibble() %>%
    dplyr::mutate(id = as.integer(id), count = as.double(count), h = as.double(h))

  # set factor variable if applicable
  if (length(name) > 0){
    for (a in 1:length(name)) {
      if (is.factor(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := factor(NA, levels = levels(covariates.comb %>% pull(sym(name[a])))))
      }
      if (is.double(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := as.double(NA))
      }
      if (is.integer(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := as.integer(NA))
      }
    }
  }

  ind <- 1
  for (j in 1:nrow(covariates.comb)) {
    ind.cov <- which(do.call(paste, covariates) == do.call(paste, covariates.comb[j, ]))
    data.sub <- Gpp$data[ind.cov, ]
    for (m in 1:Gpp$M) {
      # positions of data on line m
      ind.m <- which(data.sub$seg == m)
      y.m <- sort(as.numeric(data.sub[ind.m, ]$tp))*L$d[m]

      # bin data
      y.b <- rep(0, length(L$z[[m]]))
      for (k in 1:length(L$z[[m]])) {
        y.b[k] <- length(which(y.m < L$b[[m]][k+1] & y.m > L$b[[m]][k]))
      }

      # stack into one vector
      data$y[ind:(ind + length(y.b) - 1)] <- y.b
      data$h[ind:(ind + length(y.b) - 1)] <- L$h[m]
      ind <- ind + length(y.b)
    }
    # add bin id for every row
    data$id[((j-1)*sum(L$N.m) + 1):(j*sum(L$N.m))] <- 1:sum(L$N.m)

    # add covariates
    if (ncol(data) > 3){
      data[((j-1)*sum(L$N.m) + 1):(j*sum(L$N.m)), 4:ncol(data)] <- covariates.comb[j, ]
    }
  }
  return(data)
}

fitData <- function(Gd){
  design <- getDesign()
}

intensity <- function(X, density = FALSE){

  fit <- fitData()
}
