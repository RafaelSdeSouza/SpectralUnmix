# SpectralUnmix 0.1.0

Initial public release associated with the accepted RNAAS manuscript
`SpectralUnmix: A Torch-Based Regularized Non-negative Matrix Factorization`.

## Highlights

- torch-based regularized non-negative matrix factorization for hyperspectral
  and IFU data
- helpers for reshaping spectral cubes to spaxel-by-wavelength matrices
- standard accessors including `basis()`, `coef()`, `fitted()`, `predict()`,
  and `residuals()`
- plotting helpers for component spectra, maps, and reconstructions
- synthetic IFU demo and bundled Coelho stellar-library examples
- metadata preservation for FITS-like cube inputs and reconstructed outputs
- Quarto documentation website
