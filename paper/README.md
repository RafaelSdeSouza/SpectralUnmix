# Manuscript Drafts

This folder contains manuscript material for `SpectralUnmix`.

## Current recommendation

As of March 8, 2026, the package is better matched to a `Research Notes of the
American Astronomical Society` submission than to a full `Astronomy and
Computing` paper.

Reasons:

- the package has a coherent core API, but not yet a benchmark section
- there is not yet a public real-data vignette with astrophysical validation
- the current scientific claim is best framed as a software release / methods
  note rather than a fully evaluated software paper

## Files

- `spectralunmix_rnaas.tex`: first LaTeX draft of an RNAAS note
- `spectralunmix.bib`: bibliography for the RNAAS draft

## Suggested next steps

1. Add one real IFU example, ideally MaNGA or MUSE.
2. Add benchmark scripts comparing against classical NMF and PCA/QRPCA.
3. If those results are strong, expand this note into an `Astronomy and
   Computing` manuscript.
