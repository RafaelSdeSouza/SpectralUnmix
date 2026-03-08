#' SpectralUnmix: Spectral Unmixing for Hyperspectral and IFU Data
#'
#' Tools for reshaping spectral cubes, fitting smooth non-negative factorization
#' models, and building astronomy-focused examples for package documentation and
#' exploratory analysis.
#'
#' The package is designed around a simple workflow:
#' \enumerate{
#'   \item simulate or ingest a spectral cube;
#'   \item reshape it with [cube_to_matrix()];
#'   \item fit a low-rank unmixing model with [spectral_unmix()];
#'   \item inspect component spectra and abundance maps.
#' }
#'
#' @docType package
#' @name SpectralUnmix
NULL
