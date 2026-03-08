\name{SpectralUnmix}
\alias{SpectralUnmix}
\alias{SpectralUnmix-package}
\docType{package}
\title{SpectralUnmix: Spectral Unmixing for Hyperspectral and IFU Data}
\description{
Tools for reshaping spectral cubes, fitting smooth non-negative factorization
models, and building astronomy-focused examples for package documentation and
exploratory analysis.
}
\details{
The package is designed around a simple workflow:
\enumerate{
\item simulate or ingest a spectral cube;
\item reshape it with \code{cube_to_matrix()};
\item fit a low-rank unmixing model with \code{spectral_unmix()};
\item inspect component spectra and abundance maps.
}
}
