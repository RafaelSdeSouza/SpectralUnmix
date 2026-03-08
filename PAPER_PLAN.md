# SpectralUnmix Publication Plan

## Recommendation

As of March 8, 2026, `SpectralUnmix` is better positioned for a short research
note or software note than for a full `Astronomy and Computing` paper.

The local `qrpca_A_C/main.tex` example is closer to an `Astronomy and
Computing` submission because it already has:

- a clearly stated computational contribution
- a methodology section with algorithms
- explicit benchmark results
- a real astronomy use case on MaNGA IFU data
- package-oriented narrative throughout the manuscript

`SpectralUnmix` currently has a promising package core, but not yet the same
level of evidence or breadth.

## When To Target A&C

Aim for an `Astronomy and Computing` paper only after the package has all of
the following:

1. A stable API for fitting, plotting, reconstruction, and I/O.
2. At least one real-data astronomy case study with a public IFU cube.
3. Benchmarks against relevant baselines, not just package runtime.
4. A scientific positioning section explaining why spectral unmixing is useful
   beyond PCA/NMF boilerplate.

Without those, the paper will read like an early package announcement rather
than a software contribution paper.

## Best Immediate Path

The strongest near-term path is:

1. Finish the package into a coherent v0.1 or v0.2 release.
2. Write a short research/software note around one synthetic example and one
   public IFU example.
3. Expand that note into an `Astronomy and Computing` paper only after the
   package has benchmark and comparison sections.

## Package Work Needed

To get from the current state to a paper-ready package, prioritize:

### Core usability

- `plot_component_cube()` or `animate_component_cube()` only if channel-slice
  inspection is central to the science use case
- `select_k()` to help choose the number of components
- `reconstruct_spectrum()` for one spaxel and one component subset
- `residual_cube()` or `residual_matrix()` for model diagnostics
- `predict.spectral_unmix()` for clean reconstruction interfaces

### Astronomy workflow

- `read_ifu_fits()` to standardize FITS input and wavelength extraction
- masking support for bad pixels, NaNs, and optional variance/weight cubes
- optional continuum normalization / preprocessing helpers
- one public-data example, ideally MaNGA or MUSE

### Paper-grade evidence

- comparison against plain NMF
- comparison against PCA / QRPCA for interpretability and reconstruction
- runtime measurements for CPU and GPU
- sensitivity checks for `k`, smoothness penalty, and preprocessing choices

### Package quality

- tests for input reshaping and plotting helpers
- examples that run quickly
- one benchmark script stored in the repository
- one figure-generation script for the manuscript

## Suggested Paper Structures

### Research note / software note

Use this if you want to publish sooner.

Sections:

1. Motivation: IFU and hyperspectral decomposition problem
2. Package overview: linear mixing model and smooth NMF
3. Minimal example: synthetic IFU cube
4. Real-data illustration: one public cube
5. Limitations and future work

This version should be concise and package-centered.

### Astronomy and Computing paper

Use this only after the package is more mature.

Sections:

1. Introduction and scientific context
2. Methodology and optimization details
3. Software design and API
4. Benchmarks and scaling
5. Real IFU application
6. Comparison with PCA/NMF alternatives
7. Limitations and roadmap

## Recommended Positioning

Do not pitch `SpectralUnmix` as just "another matrix factorization package".
The stronger framing is:

"A practical spectral unmixing package for astronomy-oriented hyperspectral and
IFU analysis, focused on interpretable component spectra and abundance maps."

That is a better identity than trying to mirror PCA language too closely.

## Immediate Next Milestones

1. Add FITS-native input support.
2. Add component and residual diagnostic helpers.
3. Add one real-data vignette.
4. Add one benchmark/comparison script.
5. Draft a short note first.
