#  Copyright 2026, Rafael S. de Souza
#'
#' Spectral Unmixing for Hyperspectral and IFU Data
#'
#' Performs spectral unmixing of hyperspectral data using a non-negative matrix
#' factorization (NMF) model. The function decomposes a spectral data matrix into
#' spatial abundance maps and spectral components.
#'
#' The model assumes a linear mixing model
#'
#' \deqn{X \approx A S}
#'
#' where \eqn{X} is a matrix whose rows correspond to spectra (e.g. spaxels),
#' \eqn{A} contains spatial abundances, and \eqn{S} contains spectral
#' components. Both matrices are constrained to be non-negative.
#'
#' An optional spectral smoothness penalty may be applied to encourage
#' physically realistic spectral shapes.
#'
#' Optimization is performed using gradient-based methods implemented in
#' \code{torch}.
#'
#' @import torch
#' @export
#'
#' @param x a numeric matrix or data frame containing spectral data.
#' Each row corresponds to a spectrum (e.g. a spaxel) and columns correspond
#' to wavelength channels.
#'
#' @param k integer specifying the number of spectral components to estimate.
#'
#' @param lambda_smooth non-negative regularization parameter controlling
#' spectral smoothness. Larger values enforce smoother spectra.
#'
#' @param lr learning rate used by the optimizer.
#'
#' @param niter number of optimization iterations.
#'
#' @param center logical indicating whether column means should be subtracted.
#'
#' @param scale logical indicating whether columns should be scaled to unit
#' variance prior to factorization.
#'
#' @param cuda logical indicating whether GPU acceleration via CUDA should
#' be used.
#'
#' @return
#' \code{spectral_unmix} returns a list with class \code{"spectral_unmix"}
#' containing the following elements:
#'
#' \item{spatial}{
#' matrix of spatial abundances (or mixture weights). Each column corresponds
#' to the spatial distribution of a spectral component.
#' }
#'
#' \item{spectra}{
#' matrix of recovered spectral components. Each row represents one spectrum.
#' }
#'
#' \item{reconstruction}{
#' reconstructed data matrix obtained from the factorization \eqn{AS}.
#' }
#'
#' \item{center, scale}{
#' the centering and scaling applied to the input data, or \code{FALSE}.
#' }
#'
#' @details
#'
#' The method implements a spectral–spatial factorization of the form
#'
#' \deqn{
#' X(x,y,\lambda) =
#' \sum_{j=1}^{k} A_j(x,y) S_j(\lambda)
#' }
#'
#' which corresponds to the standard linear mixing model used in
#' hyperspectral imaging.
#'
#' In the context of integral field spectroscopy (IFU), the recovered matrices
#' can be interpreted as
#'
#' \itemize{
#' \item spectral components describing underlying stellar populations or
#' emission features;
#' \item spatial abundance maps describing the spatial distribution of each
#' component across the observed field.
#' }
#'
#' The loss function minimized during optimization is
#'
#' \deqn{
#' \|X - AS\|^2 + \lambda \|D S\|^2
#' }
#'
#' where the second term penalizes rapid spectral variations in order to
#' encourage smooth spectra.
#'
#' @references
#' Lee, D. D., & Seung, H. S. (1999). Learning the parts of objects by
#' non-negative matrix factorization. Nature, 401, 788–791.
#'
#' @examples
#' \dontrun{
#'
#' fit <- spectral_unmix(X, k = 3)
#'
#' # plot spectral components
#' matplot(t(fit$spectra), type = "l")
#'
#' # plot spatial abundance map
#' image(matrix(fit$spatial[,1], nx, ny))
#'
#' }
spectral_unmix <- function(x,
                          k = 3,
                          lambda_smooth = 0.01,
                          lr = 0.02,
                          niter = 2000,
                          center = FALSE,
                          scale = FALSE,
                          cuda = FALSE){

  if(cuda){
    device <- torch_device("cuda:0")
  } else {
    device <- torch_device("cpu")
  }

  x <- scale(x, center=center, scale=scale)

  cen <- attr(x,"scaled:center")
  sc  <- attr(x,"scaled:scale")

  X <- torch_tensor(as.matrix(x), device=device)

  p <- nrow(X)
  l <- ncol(X)

  A <- torch_rand(p,k, device=device, requires_grad=TRUE)
  S <- torch_rand(k,l, device=device, requires_grad=TRUE)

  opt <- optim_adam(list(A,S), lr=lr)

  for(i in 1:niter){

    opt$zero_grad()

    Xhat <- A$matmul(S)

    recon <- torch_mean((X - Xhat)^2)

    diff <- S[,2:l] - S[,1:(l-1)]
    smooth <- torch_mean(diff^2)

    loss <- recon + lambda_smooth*smooth

    loss$backward()
    opt$step()

    A$data()$clamp_(min=0)
    S$data()$clamp_(min=0)

  }

  A_est <- as.matrix(A$detach()$cpu())
  S_est <- as.matrix(S$detach()$cpu())

  Xhat <- A_est %*% S_est

  r <- list(
    spatial = A_est,
    spectra = S_est,
    reconstruction = Xhat,
    center = if(is.null(cen)) FALSE else cen,
    scale = if(is.null(sc)) FALSE else sc
  )

  class(r) <- "specnmf"

  return(r)
}

