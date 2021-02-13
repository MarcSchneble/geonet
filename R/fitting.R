#' Intensity Estimation on Geometric Networks with Penalized Splines
#'
#' This is the main function of the \code{geonet} package.
#' \code{intensity_pspline} estimates the intensity of a point pattern on a
#' geometric network employing penalized splines as outlined in Schneble
#' and Kauermann (2020). In distinction to
#' \code{\link[spatstat]{density.lpp}} from the \code{spatstat} package,
#' which provides kernel based tools for intensity estimation of point
#' patterns on linear networks, \code{intensity_pspline} allows to
#' incorporate covariates while also estimating the baseline intensity.
#' Covariates can be either internal or external. External covariates
#' can also be incorporated as a smooth term using penalized splines
#' with the same syntax as in \code{\link[mgcv]{gam}}.
#'
#' Details
#'
#' @param formula A one-sided formula (if a two-sided formula is supplied, the
#' left hand side of the formula is ignored). The formula can consist of either
#' linear terms as in linear models (\code{\link{lm}}) or smooth terms as
#' in \code{\link[mgcv]{gam}} formulae, where the usage is restricted to
#' smooth terms constructed with \code{\link[mgcv]{s}} and the argument
#' \code{bs} is set to
#' \code{bs = "ps"} by default, i.e. \code{intensity_pspline} can handle
#' penalized spline
#' based smooth terms. Internal linear covariates musst be supplied via
#' \code{internal()}, see examples for details.
#' @param X A point pattern on a geometric network (object of class
#' \code{gnpp}). The point pattern (specified via \code{X$data}) musst contain information on
#' all covariates included in \code{formula}.
#' @param delta The global knot distance \eqn{\delta}, a numerical vector of length one. If
#' not supplied, delta will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param h The global bin width \eqn{h}, a numerical vector of length one. If
#' not supplied, \code{h} will be chosen properly according to the geometric
#' network \code{X} which is supplied.
#' @param r The order of the penalty of the baseline intensity on the geometric
#' network, default to a penalty of order \code{r = 1}.
#' @param density \code{TRUE} if the intensity should be normalized such that it
#' can be interpreted as a density, i.e. the integral over the estimated density
#' is equal to one.
#' @return A fitted geometric network (object of class \code{gnppfit}).
#' @references Schneble, M. and G. Kauermann (2020). Intensity estimation on
#' geometric networks with penalized splines. arXiv preprint arXiv:2002.10270 .
#' @importFrom  Matrix Matrix bdiag
#' @importFrom stats update as.formula model.matrix setNames
#' @importFrom utils tail
#' @importFrom mgcv smooth.construct.ps.smooth.spec s
#' @export

intensity_pspline <- function(formula, X, delta = NULL, h = NULL, r = 1,
                             density = FALSE){
  # remove response from formula if supplied
  formula <- update(formula, NULL ~ .)

  # get linear and smooth terms
  vars <- all.vars(formula)
  if (length(vars) == 0){
    lins <- NULL
    smooths <- "G"
  } else {
    vars2 <- tail(strsplit(as.character(formula), " \\+ "),1)[[1]]
    lins <- vars[which(substring(vars2, 1, 2) != "s(")]
    smooths <- c("G", vars[which(substring(vars2, 1, 2) == "s(")])
    intern <- which(substring(vars2, 1, 9) == "internal(")
  }

  # get representation of P-splines on the network
  P <- pspline(X$network, delta, h)
  P$B <- getB(X$network, P)
  P$K <- getK(X$network, P, r)
  data <- bin_data(X, P, vars, intern)
  ind <- setNames(vector("list", length(smooths) + 1), c(smooths, "lins"))
  K <- P$K

  # design matrix of for network splines
  Z <- P$B[data$id, ]
  ind$G <- setNames(1:ncol(Z), paste0("G.", 1:ncol(Z)))

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
  fit <- fit_poisson_model(data, Z, K, ind)

  # output
  out <- list()
  out$coefficients <- setNames(fit$theta, names(unlist(ind)))
  out$V <- fit$V
  #out$effects <- effects
  out$P <- P
  out$smooth <- smooth_terms
  out$ind <- ind
  out$data <- X$data
  out$network <- X$network
  class(out) <- "gnppfit"
  out
}

#' Fit a Penalized Spline Poisson Model on a Geometric Network
#'
#' \code{fit_poisson_model} is called from \code{\link[geonet]{intensity_pspline}}
#' and performs the iterative algorithm to estimate the smoothing parameters
#' \eqn{rho} in the penalized Poisson model.
#'
#' Generalized Fellner-Schall method (Wood and Fasiolo, 2017).
#'
#' @param data soll noch weg.
#' @param Z The (sparse) model matrix where the number of coloums must
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
#' @param eps_rho The termination condition for \code{rho}.
#' @param maxit_rho Maximum number of iterations for \code{rho}.
#' @return Model fit.
#' @references Wood, S. N. and Fasiolo, M. (2017). A generalized Fellner-Schall
#'  method for smoothing parameter optimization with application to
#'  Tweedie location, scale and shape models. Biometrics 73 1071-1081.
#' @import dplyr
#' @importFrom stats optim
#' @export

fit_poisson_model <- function(data, Z, K, ind, rho = 10, rho_max = 1e5,
                     eps_rho = 0.01, maxit_rho = 100){

  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(rho, length(ind) - 1)
  Delta_rho <- Inf
  it_rho <- 0
  theta <- rep(0, ncol(Z))
  while(Delta_rho > eps_rho){
    print(rho)
    it_rho <- it_rho + 1
    theta <- scoring(theta, rho, data, Z, K, ind)

    V <- Matrix::solve(fisher(theta, rho, data, Z, K, ind))

    # updatae rho
    rho_new <- rep(NA, length(ind) - 1)
    for (i in 1:(length(ind) - 1)) {
      rho_new[i] <- as.vector(rho[i]*(Matrix::rankMatrix(K[ind[[i]], ind[[i]]], method = "qr.R")[1]/rho[i] -
                                        sum(Matrix::diag(V[ind[[i]], ind[[i]]]%*%K[ind[[i]], ind[[i]]])))/
                                (Matrix::t(theta[ind[[i]]])%*%K[ind[[i]], ind[[i]]]%*%theta[ind[[i]]]))
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
  list(theta = theta, V = V, rho = rho, it_rho = it_rho)
}
