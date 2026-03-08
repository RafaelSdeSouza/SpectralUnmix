## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SpectralUnmix)

demo <- simulate_ifu_cube()
dim(demo$cube)

## -----------------------------------------------------------------------------
x <- cube_to_matrix(demo$cube)
dim(x)

## ----eval = FALSE-------------------------------------------------------------
# fit <- spectral_unmix(
#   x,
#   k = 3,
#   lambda_smooth = 0.02,
#   niter = 300,
#   lr = 0.03
# )

## ----eval = FALSE-------------------------------------------------------------
# par(mfrow = c(1, 2))
# plot(fit, type = "spectra", wavelength = demo$wavelength)
# plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
# plot_reconstruction(fit, x, n = 4, wavelength = demo$wavelength)

