\name{coelho_demo_spectra}
\alias{coelho_demo_spectra}
\title{Demo Stellar Spectra Extracted From the Coelho Library}
\usage{
coelho_demo_spectra()
}
\value{
A list with three elements:
\itemize{
  \item \code{wavelength}: numeric wavelength grid in Angstrom.
  \item \code{matrix}: numeric matrix with one spectrum per row and one
  wavelength channel per column.
  \item \code{metadata}: data frame describing the selected spectra.
}
}
\description{
Loads four representative normalized stellar spectra extracted from the
Coelho 2014 synthetic SED library bundled in \code{extra/}. The spectra share
a common wavelength grid and are intended for lightweight demonstrations of
spectral unmixing on real library data.
}
\details{
The four bundled spectra were chosen to span distinct spectral shapes:
a cool giant, a solar-type dwarf, an A-type star, and a hot star. Each
spectrum is clipped at zero and normalized by its maximum flux so the demo
focuses on relative shape rather than absolute luminosity.
}
\examples{
demo <- coelho_demo_spectra()
dim(demo$matrix)
demo$metadata
}
