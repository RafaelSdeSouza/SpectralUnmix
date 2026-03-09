\name{coelho_stellar_subset}
\alias{coelho_stellar_subset}
\title{Curated Stellar Library Subset for NMF Demos}
\usage{
coelho_stellar_subset()
}
\value{
A list with three elements:
\itemize{
  \item \code{wavelength}: numeric wavelength grid in Angstrom.
  \item \code{matrix}: numeric matrix with one spectrum per row and one
  wavelength channel per column.
  \item \code{metadata}: data frame with spectrum identifiers, stellar type,
  and basic atmospheric parameters.
}
}
\description{
Loads a balanced subset of 100 normalized stellar spectra extracted from the
Coelho 2014 synthetic SED library. The subset contains 20 spectra from each
of five broad stellar groups and is intended for examples where a moderate
number of spectra should be represented by a small number of non-negative
components.
}
\details{
The bundled subset is balanced across five broad stellar groups:
\code{cool_giant}, \code{cool_dwarf}, \code{solar_like},
\code{a_f_star}, and \code{hot_star}. Within each group, spectra were
selected to span the available parameter range in effective temperature,
surface gravity, and metallicity. Each spectrum is clipped at zero and
normalized by its maximum flux.
}
\examples{
lib <- coelho_stellar_subset()
dim(lib$matrix)
table(lib$metadata$type)
}
