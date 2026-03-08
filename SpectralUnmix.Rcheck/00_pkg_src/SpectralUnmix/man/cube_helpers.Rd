\name{cube_to_matrix}
\alias{cube_to_matrix}
\alias{matrix_to_cube}
\alias{component_map}
\title{Helpers for Cube Reshaping and Map Extraction}
\usage{
cube_to_matrix(cube)

matrix_to_cube(x, nx, ny)

component_map(fit, nx, ny, component = 1)
}
\arguments{
\item{cube}{Numeric 3D array.}

\item{x}{Numeric matrix with one row per spatial pixel.}

\item{nx}{Number of x-axis pixels.}

\item{ny}{Number of y-axis pixels.}

\item{fit}{Result from \code{spectral_unmix()}.}

\item{component}{Component index to extract.}
}
\value{
\code{cube_to_matrix()} returns a numeric matrix with \code{nx * ny} rows and
\code{n_wave} columns.

\code{matrix_to_cube()} returns a numeric 3D array.

\code{component_map()} returns a numeric matrix of size \code{nx x ny}.
}
\description{
Utility functions for moving between cube and matrix representations and for
extracting abundance maps from a fitted model.
}
\examples{
demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
x <- cube_to_matrix(demo$cube)
rebuilt <- matrix_to_cube(x, demo$nx, demo$ny)
dim(rebuilt)
\dontrun{
fit <- spectral_unmix(x, k = 3, niter = 50)
img <- component_map(fit, demo$nx, demo$ny, 2)
}
}
