# SpectralUnmix

`SpectralUnmix` is an R package for spectral unmixing in hyperspectral imaging
and astronomical integral-field spectroscopy (IFU) cubes. It provides a compact
non-negative matrix factorization workflow with smoothness regularization,
helpers to move between cubes and matrices, and a synthetic IFU demo that can
anchor package examples, a paper, or research notes.

## Why this package

The package targets the standard linear mixing model

$$
X \approx A S
$$

where `X` contains one spectrum per spaxel, `A` stores non-negative abundance
weights, and `S` stores non-negative spectral components. In astronomy terms,
the recovered factors can be interpreted as spatial maps and component spectra
that separate continuum-like structure, emission-line dominated regions, or
other latent sources.

## Current scope

- smooth non-negative spectral unmixing implemented with `torch`
- helpers to reshape `nx x ny x n_lambda` cubes to spaxel-by-wavelength matrices
- synthetic IFU demo for reproducible examples and figures
- package website scaffold via `pkgdown`

## Installation

```r
# install.packages("remotes")
remotes::install_github("RafaelSdeSouza/SpectralUnmix")
```

`torch` must also be installed and configured in your R environment.

## Quick start

```r
library(SpectralUnmix)

demo <- simulate_ifu_cube(nx = 16, ny = 16, n_wave = 100)
fit <- spectral_unmix(demo$matrix, k = 3, niter = 200, lr = 0.03)

print(fit)

plot(fit, type = "spectra", wavelength = demo$wavelength)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
```

## Suggested workflow for a paper or research note

1. Use the synthetic IFU vignette as the reproducible baseline example.
2. Add one real-data notebook or vignette for a public cube once the target
   survey or instrument is fixed.
3. Compare `SpectralUnmix` against a classical NMF baseline and, if relevant,
   against QR/PCA-style decompositions such as the structure shown in the
   [`qrpca`](https://github.com/RafaelSdeSouza/qrpca) package.

For this repository, the strongest next step is a short methods-oriented note:
the API is still compact, so a concise research note with one synthetic example
and one real IFU application is a better fit than a full paper draft today.

## Documentation

The repository now includes:

- a getting-started vignette
- an IFU-oriented workflow vignette
- `pkgdown` configuration for a package website

Build the site locally with:

```r
pkgdown::build_site()
```

## Status

This is an early-stage research package. The core model is present, but the
next meaningful additions are benchmark examples, uncertainty diagnostics, and
one real astronomical workflow around a public IFU cube.
