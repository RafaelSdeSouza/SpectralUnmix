# SpectralUnmix

`SpectralUnmix` is an R package for smooth non-negative spectral unmixing of
hyperspectral data and astronomical integral-field spectroscopy (IFU) cubes.
The package fits a low-rank model that decomposes a spectral matrix into
non-negative component spectra and spatial abundance maps.

The model is

$$
X \approx A S
$$

where `X` is a spaxel-by-wavelength matrix, `A` contains spatial abundances,
and `S` contains component spectra.

## Installation

```r
# install.packages("remotes")
remotes::install_github("RafaelSdeSouza/SpectralUnmix")
```

`torch` must be installed and configured in the local R environment.

## Basic usage

```r
library(SpectralUnmix)

demo <- simulate_ifu_cube(nx = 16, ny = 16, n_wave = 100)

fit <- spectral_unmix(
  demo$matrix,
  k = 3,
  lambda_smooth = 0.02,
  niter = 300,
  lr = 0.03
)

print(fit)
summary(fit)
```

## Accessors

The fitted object supports a compact interface inspired by `prcomp` and `NMF`.

```r
# abundance matrix (spaxels x components)
fit$spatial

# component spectra (components x wavelength)
fit$spectra

# NMF-style accessors
basis(fit)        # wavelength x components
coef(fit)         # components x spaxels

# fitted values and residuals
fitted(fit)
residuals(fit, x = demo$matrix)

# in-sample or new-data prediction
predict(fit)
predict(fit, newdata = demo$cube, type = "spatial")

# metadata carried by the fit
cube_metadata(fit)
```

## Visualization

```r
plot(fit, type = "spectra", wavelength = demo$wavelength)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
plot(fit, type = "loss")
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
```

## IFU workflow

```r
library(FITSio)

X <- readFITS("cube.fits")
Mat <- cube_to_matrix(X)
Mat[!is.finite(Mat)] <- 0
Mat <- pmax(Mat, 0)

fit <- spectral_unmix(Mat, k = 5, lr = 0.01, niter = 5000)

maps <- predict(fit, type = "spatial")
cube_hat <- predict(
  fit,
  type = "cube",
  nx = dim(X$imDat)[1],
  ny = dim(X$imDat)[2]
)

# recover stored FITS-side metadata such as headers or redshift if present
meta <- cube_metadata(cube_hat)
```

When `cube_to_matrix()` receives a FITS-like list object, it now carries
non-image entries such as headers and other metadata through the matrix, the
fitted object, and reconstructed cubes.

## Real spectra demo

```r
real_demo <- coelho_demo_spectra()
dim(real_demo$matrix)

fit_real <- spectral_unmix(
  real_demo$matrix,
  k = 3,
  lambda_smooth = 0.001,
  niter = 400,
  lr = 0.03
)
```

## Stellar library subset

```r
stellar_lib <- coelho_stellar_subset()
dim(stellar_lib$matrix)
table(stellar_lib$metadata$type)

fit_lib <- spectral_unmix(
  stellar_lib$matrix,
  k = 5,
  lambda_smooth = 0.001,
  niter = 500,
  lr = 0.03
)
```

## Documentation

Package articles are provided in the `vignettes/` directory.

Website: [https://rafaelsdesouza.github.io/SpectralUnmix/](https://rafaelsdesouza.github.io/SpectralUnmix/)
