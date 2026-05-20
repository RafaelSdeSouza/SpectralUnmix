#' Demo Stellar Spectra Extracted From the Coelho Library
#'
#' Loads four representative normalized stellar spectra extracted from the
#' Coelho 2014 synthetic SED library bundled in `extra/`. The spectra share a
#' common wavelength grid and are intended for lightweight demonstrations of
#' spectral unmixing on real library data.
#'
#' @return A list with three elements:
#' \itemize{
#'   \item `wavelength`: numeric wavelength grid in Angstrom.
#'   \item `matrix`: numeric matrix with one spectrum per row and one
#'     wavelength channel per column.
#'   \item `metadata`: data frame describing the selected spectra.
#' }
#'
#' @details
#' The four bundled spectra were chosen to span distinct spectral shapes:
#' a cool giant, a solar-type dwarf, an A-type star, and a hot star. Each
#' spectrum is clipped at zero and normalized by its maximum flux so the demo
#' focuses on relative shape rather than absolute luminosity.
#'
#' @examples
#' demo <- coelho_demo_spectra()
#' dim(demo$matrix)
#' demo$metadata
#' @export
coelho_demo_spectra <- function() {
  spectra_path <- system.file("extdata", "coelho-demo-spectra.csv", package = "SpectralUnmix")
  metadata_path <- system.file("extdata", "coelho-demo-metadata.csv", package = "SpectralUnmix")

  if (!nzchar(spectra_path) || !nzchar(metadata_path)) {
    stop("Bundled demo spectra were not found in inst/extdata.", call. = FALSE)
  }

  spectra_df <- utils::read.csv(spectra_path, check.names = FALSE)
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)

  wavelength <- as.numeric(spectra_df[[1L]])
  spectra <- t(as.matrix(spectra_df[, -1L, drop = FALSE]))
  rownames(spectra) <- names(spectra_df)[-1L]
  colnames(spectra) <- wavelength

  list(
    wavelength = wavelength,
    matrix = spectra,
    metadata = metadata
  )
}

#' Curated Stellar Library Subset for NMF Demos
#'
#' Loads a balanced subset of 100 normalized stellar spectra extracted from the
#' Coelho 2014 synthetic SED library. The subset contains 20 spectra from each
#' of five broad stellar groups and is intended for examples where a moderate
#' number of spectra should be represented by a small number of non-negative
#' components.
#'
#' @return A list with three elements:
#' \itemize{
#'   \item `wavelength`: numeric wavelength grid in Angstrom.
#'   \item `matrix`: numeric matrix with one spectrum per row and one
#'     wavelength channel per column.
#'   \item `metadata`: data frame with spectrum identifiers, stellar type, and
#'     basic atmospheric parameters.
#' }
#'
#' @details
#' The bundled subset is balanced across five broad stellar groups:
#' `cool_giant`, `cool_dwarf`, `solar_like`, `a_f_star`, and `hot_star`.
#' Within each group, spectra were selected to span the available parameter
#' range in effective temperature, surface gravity, and metallicity. Each
#' spectrum is clipped at zero and normalized by its maximum flux.
#'
#' @examples
#' lib <- coelho_stellar_subset()
#' dim(lib$matrix)
#' table(lib$metadata$type)
#' @export
coelho_stellar_subset <- function() {
  subset_path <- system.file("extdata", "coelho-stellar-subset.rds", package = "SpectralUnmix")

  if (!nzchar(subset_path)) {
    stop("Bundled stellar subset was not found in inst/extdata.", call. = FALSE)
  }

  readRDS(subset_path)
}
