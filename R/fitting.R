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
#' @param formula A one-sided formula (if a two-sided formula is supplied, the
#' left hand side of the formula is ignored). The formula can consist of either
#' linear terms as in linear models (\code{\link{lm}}) or smooth terms as
#' in \code{\link[mgcv]{gam}} formulae, where the usage is restricted to
#' smooth terms constructed with \code{\link[mgcv]{s}} and the argument
#' \code{bs} is set to
#' \code{bs = "ps"} by default, i.e. \code{intensity_pspline} can handle
#' penalized spline
#' based smooth terms. Internal linear covariates must be supplied via
#' \code{internal()}, see examples for details.
#' @param delta The global knot distance \eqn{\delta}, a numerical vector of length one. If
#' not supplied, delta will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param h The global bin width \eqn{h}, a numerical vector of length one. If
#' not supplied, \code{h} will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param r The order of the penalty of the baseline intensity on the geometric
#' network, default to a penalty of order \code{r = 1}.
#' @param scale A named list which specifies the rescaling of network related
#' covariates. Currently, only x- and y-coordinates can be scaled.
#' @param density \code{TRUE} if the intensity should be normalized such that it
#' can be interpreted as a density, i.e. the integral over the estimated density
#' is equal to one.
#' @param verbose If \code{TRUE}, prints information on the process of the fitting
#' algorithm.
#' @return A fitted geometric network (object of class \code{gnppfit}).
#' @references Schneble, M. and G. Kauermann (2020). Intensity estimation on
#' geometric networks with penalized splines. arXiv preprint arXiv:2002.10270 .
#' @importFrom  Matrix Matrix bdiag
#' @importFrom stats update as.formula model.matrix setNames
#' @importFrom utils tail
#' @importFrom mgcv smooth.construct.ps.smooth.spec s
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export
#' @examples
#' library(geonet)
#' X <- runifgn(50, small_gn)
#' delta <- 0.2
#' h <- 0.1
#' r <- 2
#' model <- intensity_pspline(X, delta = delta, h = h, r = r)
#' summary(model)
#' plot(model)

intensity_pspline <- function(X, formula = ~1, delta = NULL, h = NULL, r = 1,
                             scale = NULL, density = FALSE, verbose = FALSE){
  # remove response from formula if supplied
  formula <- update(formula, NULL ~ .)

  # get linear and smooth terms
  vars <- all.vars(formula)
  if (length(vars) == 0){
    lins <- NULL
    smooths <- "G"
    intern <- NULL
  } else {
    vars2 <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
    lins <- vars[which(substring(vars2, 1, 2) != "s(")]
    smooths <- c("G", vars[which(substring(vars2, 1, 2) == "s(")])
    intern <- which(substring(vars2, 1, 9) == "internal(")
  }

  # get representation of P-splines on the network
  knots <- network_knots(X$network, delta)
  bins <- network_bins(X$network, h)
  #P <- list(splines = knots, bins = bins)
  data <- bin_data(X, bins, vars, intern, scale)
  ind <- setNames(vector("list", length(smooths) + 1), c(smooths, "lins"))

  # design matrix of for network splines
  B <- bspline_design(X$network, knots, bins)
  Z <- B[data$id, ]
  ind$G <- setNames(1:ncol(Z), paste0("G.", 1:ncol(Z)))
  K <- network_penalty(X$network, knots, r)

  smooth_terms <- setNames(vector("list", length(smooths) - 1), smooths[-1])
  # design for smooth terms
  if (length(smooths) > 1) {
    for (i in 2:length(smooths)) {
      sm <- smooth.construct.ps.smooth.spec(eval(parse(text = vars2[i - 1])),
                                            data = data, knots = NULL)
      D <- sm$D[, -1]
      K <- bdiag(K, t(D)%*%D)
      Z_smooth <- sm$X
      smooth_terms[[smooths[i]]]$range <- range(data[[smooths[i]]])
      smooth_terms[[smooths[i]]]$knots <- sm$knots
      smooth_terms[[smooths[i]]]$l <- sm$m[1] + 1
      smooth_terms[[smooths[i]]]$r <- sm$m[2]
      smooth_terms[[smooths[i]]]$ident <- colMeans(Z_smooth)
      Z <- cbind(Z, sweep(Z_smooth, 2, smooth_terms[[smooths[i]]]$ident)[, -1])
      ind[[smooths[i]]] <- setNames((ncol(Z) - ncol(Z_smooth) + 2):ncol(Z),
                                    paste0(smooths[i], ".", 1:(ncol(Z_smooth) - 1)))
    }
  }

  # design for linear terms
  if (length(lins) > 0){
    formula_lins <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model_matrix_lins <- model.matrix(formula_lins, data = data)[, -1]
    ind[["lins"]] <- setNames((ncol(Z) + 1):(ncol(Z) + ncol(model_matrix_lins)),
                              colnames(model_matrix_lins))
    Z <- cbind(Z, model_matrix_lins)
    K <- bdiag(K, Matrix(0, ncol(model_matrix_lins), ncol(model_matrix_lins), sparse = TRUE))
  }

  # fit the model
  if (verbose) cat("Finished preprocessing. Start fitting the model.\n")
  fit <- fit_poisson_model(data, Z, K, ind, verbose = verbose)

  # output
  out <- list()
  out$coefficients <- setNames(fit$theta, names(unlist(ind)))
  out$V <- fit$V
  #out$effects <- effects
  #out$P <- P
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
#' @param rho An initial estimate of the smoothing parameter. Either a vector
#' of length one or a vector which corresponds to the count of smooth terms
#' (including the baseline intensity) in the model.
#' @param rho_max If \code{rho} exceeds \code{rho_max}, the algorithm stops
#' and returns the current values.
#' @param eps_theta The termination condition.
#' @param maxit_rho Maximum number of iterations for \code{rho}.
#' @param verbose If \code{TRUE}, prints information on the process of the fitting
#' algorithm.
#' @return Model fit.
#' @references Wood, S. N. and Fasiolo, M. (2017). A generalized Fellner-Schall
#'  method for smoothing parameter optimization with application to
#'  Tweedie location, scale and shape models. Biometrics 73 1071-1081.
#' @import dplyr
#' @importFrom stats optim
#' @author Marc Schneble \email{marc.schneble@@stat.uni-muenchen.de}
#' @export

fit_poisson_model <- function(data, Z, K, ind, rho = 10, rho_max = 1e5,
                     eps_theta = 0.001, maxit_rho = 100, verbose = FALSE){

  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(rho, length(ind) - 1)
  Delta_theta <- Inf
  it_rho <- 0
  theta <- rep(0, ncol(Z))
  duration <- NULL
  while(tail(Delta_theta, 1) > eps_theta){
    print(rho)
    if (verbose) start <- Sys.time()
    it_rho <- it_rho + 1
    if (it_rho == 1) {
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

    if (any(rho_new > rho_max)) break
    if (any(rho_new < 0)){
      warning("rho = 0 has occurred")
    }
    if (it_rho > maxit_rho){
      warning("Stopped estimation of rho because maximum number of iterations has been reached!")
      break
    }
    Delta_theta <- c(Delta_theta, sqrt(sum((theta_new - theta)^2))/sqrt(sum((theta)^2)))
    if (verbose & it_rho >= 5) {
      duration <- c(duration, as.vector(difftime(Sys.time(), start, units = "mins")))
      Delta_theta_recent <- tail(Delta_theta, 4)
      mean_reduction <- mean((abs(Delta_theta_recent[-1] -
                                    Delta_theta_recent[-4])/
                                Delta_theta_recent[-4]))
      k <- ceiling(log(eps_theta/tail(Delta_theta, 1))/log(1 - mean_reduction))
      cat(paste("Expected number of further iterations:", k,
                "(expected duration:", round(mean(tail(duration, 4))*k, 1), "minutes)", "\n"))
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

  list(theta = theta, V = V, rho = rho, it_rho = it_rho, edf = edf)
}
