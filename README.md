# SpectralUnmix

**SpectralUnmix** is an R package for spectral–spatial decomposition of hyperspectral datasets and integral field spectroscopy (IFU) cubes.

The package implements spectral unmixing methods that decompose a spectral cube into a small number of physically interpretable spectral components and spatial abundance maps.

---

## Model

SpectralUnmix assumes a linear mixing model

$$
X(x,y,\lambda) = \sum_{j=1}^{k} A_j(x,y) S_j(\lambda)
$$

where

- $$S_j(\lambda)$$ are **spectral components**
- $$A_j(x,y)$$ are **spatial abundance maps**

In matrix form
$X \approx AS$
where

| Matrix | Meaning |
|------|------|
| **X** | spectra per spaxel |
| **A** | spatial abundance maps |
| **S** | spectral components |

---

## Features

- spectral unmixing for hyperspectral and IFU data
- non-negative matrix factorization
- optional spectral smoothness regularization
- GPU acceleration via **torch**
- scalable to large spectral datasets

---

## Installation

```r
devtools::install_github("rsdesouza/SpectralUnmix")


library(SpectralUnmix)

fit <- spectral_unmix(X, k = 3)

# spectral components
matplot(t(fit$spectra), type = "l")

# spatial abundance maps
image(matrix(fit$spatial[,1], nx, ny))

