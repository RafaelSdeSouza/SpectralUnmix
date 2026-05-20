# Zenodo Release Checklist for SpectralUnmix

This checklist is tailored to the accepted manuscript:

**SpectralUnmix: A Torch-Based Regularized Non-negative Matrix Factorization**

## 1. Freeze the release version

Pick the exact package version to cite in the paper, then make sure
`DESCRIPTION`, `NEWS.md` if used, and the GitHub release tag all match.

Example:

```bash
git tag -a v0.1.0 -m "SpectralUnmix v0.1.0"
git push origin v0.1.0
```

If you need a paper-specific release, create that tag only after the accepted
manuscript fixes are complete.

## 2. Check release metadata before Zenodo import

The repository now includes `.zenodo.json`, which Zenodo will read.

Before publishing the Zenodo record, confirm:

- title is correct
- creator names are real names, not GitHub handles
- ORCIDs are present
- affiliations are correct
- license matches GitHub and the package license

## 3. Connect GitHub to Zenodo

In Zenodo:

1. Log in
2. Open the GitHub integration page
3. Enable the `SpectralUnmix` repository

Once enabled, Zenodo will archive each tagged GitHub release.

## 4. Create the GitHub release

On GitHub:

1. Open `Releases`
2. Create a new release from the frozen tag
3. Use the same version number as the package and paper citation

Zenodo should then generate:

- a **version DOI** for that exact frozen release
- a **concept DOI** for the repository as a whole

For the paper and the editor, send the **version DOI**.

## 5. Submit the Zenodo record to the AAS community

After Zenodo creates the deposit:

1. Submit it to the AAS Journals community
2. Include a message to the curators with:
   - manuscript title
   - full author name
   - eJP submission ID: `#AAS74884`

## 6. Update the paper and README

Once the version DOI exists:

- add the DOI citation to the final article text and bibliography
- add the DOI badge or citation line to the package README
- keep the GitHub release tag, Zenodo record, and package version aligned

## 7. Notebook check

If notebooks are added to the repository before final publication, review the
AAS notebook guidance before linking them in the paper.
