## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SpectralUnmix)

demo <- simulate_ifu_cube(nx = 14, ny = 14, n_wave = 90)
dim(demo$matrix)

## ----eval = FALSE-------------------------------------------------------------
# fit <- spectral_unmix(demo$matrix, k = 3, niter = 250, lr = 0.03)
# print(fit)
# plot(fit, type = "spectra", wavelength = demo$wavelength)
# plot(fit, type = "loss")

## ----eval = FALSE-------------------------------------------------------------
# map1 <- component_map(fit, nx = demo$nx, ny = demo$ny, component = 1)
# image(map1, main = "Component 1", xlab = "x", ylab = "y")
# plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)

