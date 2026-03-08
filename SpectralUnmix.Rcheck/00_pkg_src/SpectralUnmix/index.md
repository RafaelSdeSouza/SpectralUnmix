---
title: SpectralUnmix
---

`SpectralUnmix` is built for low-rank decomposition of hyperspectral cubes and
astronomical IFU data. The package centers the workflow around one practical
question: can a spectral cube be decomposed into a small number of smooth,
non-negative spectral components and spatial abundance maps that remain easy to
interpret?

The current package is intentionally compact. It provides:

- a smooth NMF-style model implemented with `torch`
- helpers for reshaping between cubes and spaxel matrices
- a synthetic IFU generator for examples, demos, and figure prototyping

This structure is suitable for a research note or a methods-focused package
paper once a real-data case study is added.
