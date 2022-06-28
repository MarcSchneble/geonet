#' Intensity Estimation on Geometric Networks with Penalized Splines
#'
#' This is the main function of the \code{geonet} package.
#' \code{intensity_pspline} estimates the intensity of a point pattern on a
#' geometric network employing penalized splines as outlined in Schneble
#' and Kauermann (2020). In distinction to
#' \code{\link[spatstat.linnet]{density.lpp}} from the \code{spatstat.linnet} package,
#' which provides kernel based tools for intensity estimation of point
#' patterns on linear networks, \code{intensity_pspline} allows to
#' incorporate covariates while also estimating the baseline intensity.
#' Covariates can be either internal or external. External covariates
#' can also be incorporated as a smooth term using penalized splines
#' with the same syntax as in \code{\link[mgcv]{gam}}.
#'
#' @param X A point pattern on a geometric network (object of class
#' \code{gnpp}). The data (\code{X$data}) must contain information on
#' all covariates included in \code{formula}.
#' @param ... Other arguments. The following arguments must match exactly.
#' @param formula A one-sided formula (if a two-sided formula is supplied, the
#' left hand side of the formula is ignored). The formula can consist of either
#' linear terms as in linear models (\code{\link{lm}}) or smooth terms as
#' in \code{\link[mgcv]{gam}} formulae, where the usage is restricted to
#' smooth terms constructed with \code{\link[mgcv]{s}} and the argument
#' \code{bs} is set to
#' \code{bs = "ps"} by default, i.e. \code{intensity_pspline} can handle
#' penalized spline
#' based smooth terms.
#' @param delta The global knot distance \eqn{\delta}, a numerical vector of length one. If
#' not supplied, delta will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param h The global bin width \eqn{h}, a numerical vector of length one. If
#' not supplied, \code{h} will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param r The order of the penalty of the baseline intensity on the geometric
#' network, default to a penalty of order \code{r = 2}.
#' @param scale A named list which specifies the rescaling of network related
#' covariates. Currently, only x- and y-coordinates can be scaled.
#' @param density \code{TRUE} if the intensity should be normalized such that it
#' can be interpreted as a density, i.e. the integral over the estimated density
#' is equal to one.
#' @param verbose If \code{TRUE}, prints information on the process of the fitting
#' algorithm.
#' @param control A list of optional arguments which control the convergence
#' of the fitting algorithm. See "Details".
#' @return A fitted geometric network (object of class \code{gnppfit}).
#' @references Schneble, M. and G. Kauermann (2020). Intensity estimation on
#' geometric networks with penalized splines. arXiv preprint arXiv:2002.10270 .
#' @importFrom Matrix Matrix bdiag
#' @importFrom stats update as.formula model.matrix setNames
#' @importFrom utils tail
#' @importFrom mgcv smooth.construct.ps.smooth.spec s
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export
#' @examples
#' library(geonet)
#' X <- runifgn(50, small_gn)
#' delta <- 0.2
#' model <- intensity_pspline(X, delta = delta)
#' summary(model)
#' plot(model)

intensity_pspline <- function(X, ..., formula = ~1, delta = "0", h = "0.5", r = 2,
                              scale = NULL, density = FALSE, verbose = FALSE,
                              control = list()){

  if (density) stop("Computing the density is currently not supported.")
  stopifnot(class(X) == "gnpp" & is.character(delta) & is.character(h) &
              is.logical(verbose) & r %in% c(1, 2))

  # remove response from formula if supplied
  formula <- update(formula, NULL ~ .)

  # get linear and smooth terms
  vars <- all.vars(formula)
  if (length(vars) == 0){
    vars_lins <- NULL
    vars_smooths <- "G"
    vars_internal <- NULL
    ind_smooths <- NULL
  } else {
    vars2 <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
    if (length(setdiff(vars, c("x", "y", "dist2V",
                               names(summary(X)$covariates$internal),
                               names(summary(X)$covariates$external)))) > 0) {
      stop("At least one covariate was not found in the data!")
    }
    formula_string <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
    vars_internal <- intersect(vars, c("x", "y", "dist2V",
                                       names(summary(X)$covariates$internal)))
    vars_lins <- vars[which(substring(formula_string, 1, 2) != "s(")]
    ind_smooths <- which(substring(formula_string, 1, 2) == "s(")
    vars_smooths <- c("G", vars[ind_smooths])
    if (length(intersect(vars_internal, vars_smooths)) > 0) {
      stop("Internal covariates can only be included as linear terms!")
    }
  }

  # get representation of P-splines on the network
  setup <- delta_h_global(X$network, delta, h)
  knots <- network_knots(X$network, setup$delta)
  bins <- network_bins(X$network, setup$h)

  data <- bin_data(X, bins = bins, vars = vars, vars_internal = vars_internal,
                   scale = scale)
  ind <- setNames(vector("list", length(vars_smooths) + 1),
                  c(vars_smooths, "lins"))

  # design matrix of for network splines
  B <- bspline_design(X$network, knots, bins)
  Z <- B[data$id, ]
  ind$G <- setNames(1:ncol(Z), paste0("G.", 1:ncol(Z)))
  K <- network_penalty(X$network, knots, r)

  smooth_terms <- setNames(vector("list", length(vars_smooths) - 1),
                           vars_smooths[-1])
  # design for smooth terms
  if (length(ind_smooths) > 0) {
    for (i in 1:length(ind_smooths)) {
      sm <- smooth.construct.ps.smooth.spec(
        eval(parse(text = formula_string[ind_smooths][i])),
                                            data = data, knots = NULL)
      D <- sm$D[, -1]
      K <- bdiag(K, t(D)%*%D)
      Z_smooth <- sm$X
      var_name <- vars[ind_smooths][i]
      smooth_terms[[var_name]]$range <- range(data[[var_name]])
      smooth_terms[[var_name]]$knots <- sm$knots
      smooth_terms[[var_name]]$l <- sm$m[1] + 1
      smooth_terms[[var_name]]$r <- sm$m[2]
      smooth_terms[[var_name]]$ident <- colMeans(Z_smooth)
      Z <- cbind(Z, sweep(Z_smooth, 2, smooth_terms[[var_name]]$ident)[, -1])
      ind[[var_name]] <- setNames((ncol(Z) - ncol(Z_smooth) + 2):ncol(Z),
                                    paste0(var_name, ".", 1:(ncol(Z_smooth) - 1)))
    }
  }

  # design for linear terms
  if (length(vars_lins) > 0){
    formula_lins <- as.formula(paste("~", paste(vars_lins, collapse = " + ")))
    model_matrix_lins <- as.matrix(model.matrix(formula_lins, data = data)[, -1])
    ind[["lins"]] <- setNames((ncol(Z) + 1):(ncol(Z) + ncol(model_matrix_lins)),
                              colnames(model_matrix_lins))
    Z <- cbind(Z, model_matrix_lins)
    K <- bdiag(K, Matrix(0, ncol(model_matrix_lins), ncol(model_matrix_lins), sparse = TRUE))
  }

  # fit the model
  if (verbose) cat("Finished preprocessing. Start fitting the model.\n")
  fit <- fit_poisson_model(data, Z, K, ind, verbose = verbose,
                           control = control)

  # output
  out <- list()
  out$coefficients <- setNames(fit$theta, names(unlist(ind)))
  out$V <- fit$V
  #out$effects <- effects
  #out$P <- P
  out$setup <- setup
  out$knots <- knots
  out$bins <- bins
  out$smooth <- smooth_terms
  out$ind <- ind
  out$data <- X$data
  out$network <- X$network
  out$formula <- formula
  out$it_rho <- fit$it_rho
  out$edf <- fit$edf
  class(out) <- "gnppfit"
  out
}

#' Fit a Penalized Spline Poisson Model on a Geometric Network
#'
#' \code{fit_poisson_model} is called from \code{\link[geonet]{intensity_pspline}}
#' and performs the iterative algorithm to estimate the model parameters and the
#'  smoothing parameters \eqn{rho} in the penalized Poisson model.
#'
#' Smoothing parameters are estimated using the generalized Fellner-Schall
#' method (Wood and Fasiolo, 2017).
#'
#' @param data The binned data.
#' @param Z The (sparse) model matrix where the number of columns must
#' correspond to the length of the vector of model coefficients \code{theta}.
#' @param K A (sparse) square penalty matrix of with the same dimension as
#' \code{theta}.
#' @param ind A list which contains the indices belonging to each smooth term
#' and the linear terms.
#' @param verbose If \code{TRUE}, prints information on the process of the fitting
#' algorithm.
#' @param control A list of optional arguments which control the convergence
#' of the fitting algorithm. See "Details".
#' @return Model fit.
#' @references Wood, S. N. and Fasiolo, M. (2017). A generalized Fellner-Schall
#'  method for smoothing parameter optimization with application to
#'  Tweedie location, scale and shape models. Biometrics 73 1071-1081.
#' @import dplyr
#' @importFrom stats optim
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}

fit_poisson_model <- function(data, Z, K, ind, verbose = FALSE,
                              control = list()){

  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(ifelse(!is.null(control$rho_start), control$rho_start, 10),
             length(ind) - 1)
  rho_max <- ifelse(!is.null(control$rho_max), control$rho_max, 1e5)
  it_max <- ifelse(!is.null(control$it_max), control$it_max, 100)
  eps_theta <- ifelse(!is.null(control$eps_theta), control$eps_theta, 1e-3)

  Delta_theta <- Inf
  it <- 0
  theta <- rep(0, ncol(Z))
  duration <- NULL
  while(tail(Delta_theta, 1) > eps_theta){
    if (verbose) start <- Sys.time()
    it <- it + 1
    if (it == 1) {
      theta_new <- scoring(theta, rho, data, Z, K, ind)
    } else {
      theta_new <-
        as.vector(theta + Matrix::solve(fisher(theta, rho, data, Z, K, ind))%*%
                    score(theta, rho, data, Z, K, ind))

    }
    V <- Matrix::solve(fisher(theta_new, rho, data, Z, K, ind))

    # update rho
    rho_new <- rep(NA, length(ind) - 1)
    for (i in 1:(length(ind) - 1)) {
      rho_new[i] <- as.vector(rho[i]*(Matrix::rankMatrix(K[ind[[i]], ind[[i]]], method = "qr.R")[1]/rho[i] -
                                        sum(Matrix::diag(V[ind[[i]], ind[[i]]]%*%K[ind[[i]], ind[[i]]])))/
                                (Matrix::t(theta_new[ind[[i]]])%*%K[ind[[i]], ind[[i]]]%*%theta_new[ind[[i]]]))
    }

    if (any(rho_new > rho_max)) {
      warning("Stopped estimation of rho because maximum rho has been reached!")
      break
    }
    if (any(rho_new < 0)){
      warning("rho = 0 has occurred")
    }
    if (it > it_max){
      warning("Stopped estimation of rho because maximum number of iterations has been reached!")
      break
    }
    Delta_theta <- c(Delta_theta, sqrt(sum((theta_new - theta)^2))/sqrt(sum((theta)^2)))
    if (verbose & it >= 5) {
      duration <- c(duration, as.vector(difftime(Sys.time(), start, units = "mins")))
      Delta_theta_recent <- tail(Delta_theta, 4)
      mean_reduction <- mean((abs(Delta_theta_recent[-1] -
                                    Delta_theta_recent[-4])/
                                Delta_theta_recent[-4]))
      k <- ceiling(log(eps_theta/tail(Delta_theta, 1))/log(1 - mean_reduction))
      cat(paste("Expected number of further iterations:", k,
                "(expected duration:", round(mean(tail(duration, 4))*k, 2), "minutes)", "\n"))
    }
    rho <- rho_new
    theta <- theta_new
  }
  # degrees of freedom
  mu <- exp(as.vector(Z%*%theta) + log(data$h) + log(data$offset))
  Mu <- Matrix::Matrix(0, nrow(Z), nrow(Z))
  diag(Mu) <- mu
  H <- Matrix::solve(fisher(theta, rho, data, Z, K, ind))%*%(Matrix::t(Z)%*%Mu%*%Z)
  edf <- Matrix::diag(H)

  list(theta = theta, V = V, rho = rho, it_rho = it, edf = edf)
}



#' Intensity Estimation on Geometric Networks based on Kernel Smoothing
#'
#' @param X A point pattern on a geometric network (object of class
#' \code{gnpp}).
#' @param kernel If \code{kernel = "heat"}, a heat kernel is used. If
#' \code{kernel = "Euclidean"}, a two-dimensional kernel smoother is used.
#'
#' @return A fitted point process on a linear network, an object of class
#' \code{lppfit}.
#' @import spatstat.linnet
#' @import spatstat.core
#' @export
#' @examples
#' X <- runifgn(n = 50, G = small_gn)
#' fit <- intensity_kernel(X)
#' plot(fit)


intensity_kernel <- function(X, kernel = "heat") {
  stopifnot(class(X) == "gnpp")
  match.arg(arg = kernel, choices = c("heat", "Euclidean"))
  Y <- as_lpp(X)
  if (kernel == "heat") {
    sigma <- bw.lppl(Y, distance = "path")
    fit <- density.lpp(Y, sigma = as.numeric(sigma))
  }
  if (kernel == "Euclidean") {
    sigma <- bw.scott.iso(Y)
    fit <- density.lpp(Y, sigma = sigma, distance = "euclidean")
  }
  fit$network <- X$network
  class(fit) <- c("lppfit", class(fit))
  fit
}
