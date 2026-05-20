# Release Review for RNAAS / Zenodo

This file is the final pre-release review checklist for the frozen software
version linked to manuscript `#AAS74884`.

## Release target

- package version: `0.1.0`
- proposed tag: `v0.1.0`
- manuscript title:
  `SpectralUnmix: A Torch-Based Regularized Non-negative Matrix Factorization`

## Files to review before tagging

- [DESCRIPTION](DESCRIPTION)
- [.zenodo.json](.zenodo.json)
- [CITATION.cff](CITATION.cff)
- [LICENSE](LICENSE)
- [README.md](README.md)
- [NEWS.md](NEWS.md)
- [ZENODO_RELEASE.md](ZENODO_RELEASE.md)

## Tag and release commands

Run only after confirming the repository is in its final accepted state:

```bash
git status
git tag -a v0.1.0 -m "SpectralUnmix v0.1.0"
git push origin v0.1.0
```

Then create the GitHub Release from tag `v0.1.0`.

## Zenodo review points

- creator names use real names
- ORCID is present
- affiliation is correct
- title matches manuscript
- license matches GitHub repository
- use the **version DOI** in the paper and in the reply to the data editor

## Suggested GitHub release title

`SpectralUnmix v0.1.0`

## Suggested GitHub release notes

Initial frozen release accompanying the accepted RNAAS manuscript
`SpectralUnmix: A Torch-Based Regularized Non-negative Matrix Factorization`.

This release includes the torch-based regularized NMF implementation, cube
reshaping helpers, plotting/accessor methods, metadata-preserving IFU outputs,
bundled demo datasets, and the Quarto documentation website.
