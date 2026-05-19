# SpectralUnmix v2 Design Note

This note describes a more general modeling direction for `SpectralUnmix`
without changing the current package behavior yet.

The central idea is to treat spectral unmixing as a latent-variable model with:

- non-negative latent factors
- signed observed data
- an explicit observation model
- optional regularization on the latent spectra and abundances

This is broader than classical NMF and cleaner than building the package around
one special-purpose update rule.

## Motivation

The current package already fits

\[
X \approx AS
\]

with:

- `A >= 0`: abundance or weight matrix
- `S >= 0`: component spectra

while allowing the observed matrix `X` to contain negative values after
preprocessing. That is already useful for astronomy, where the true signal is
non-negative but the measured data may be negative because of noise,
background subtraction, or calibration residuals.

What is still missing is a more explicit statistical framing of that problem.

## Proposed model

Let:

- `X ∈ R^{n × p}` be the observed matrix
- `A ∈ R_+^{n × k}` be the latent non-negative abundance matrix
- `S ∈ R_+^{k × p}` be the latent non-negative component spectra

The observation model is:

\[
X = AS + \varepsilon
\]

where `\varepsilon` is allowed to be signed.

This gives the right separation of roles:

- the latent signal is constrained to be non-negative
- the observations are allowed to fluctuate above and below that signal

## Torch parameterization

In a torch-based implementation, positivity should ideally be built into the
parameterization rather than enforced only by post-update clamping.

Introduce unconstrained latent variables:

\[
U_A \in \mathbb{R}^{n \times k}, \quad U_S \in \mathbb{R}^{k \times p}
\]

and map them to positive factors using:

\[
A = \operatorname{softplus}(U_A) + \epsilon
\]

\[
S = \operatorname{softplus}(U_S) + \epsilon
\]

for a small `\epsilon > 0`.

This has a few advantages over hard clamping:

- positivity is intrinsic to the model
- gradients stay smooth
- extensions to weighted and robust likelihoods are cleaner

## Observation models

The main generalization should come from the likelihood or loss, not from a
single solver-specific update rule.

### 1. Gaussian

Default when no weights are available:

\[
X_{ij} \sim \mathcal{N}((AS)_{ij}, \sigma^2)
\]

Loss:

\[
\mathcal{L}_{\mathrm{gaussian}} =
\frac{1}{2}\sum_{ij}(X_{ij} - (AS)_{ij})^2
\]

### 2. Weighted Gaussian

Recommended astronomy-facing default when inverse-variance information is
available:

\[
X_{ij} \sim \mathcal{N}((AS)_{ij}, \sigma_{ij}^2)
\]

with weights `w_{ij} = \sigma_{ij}^{-2}` and loss

\[
\mathcal{L}_{\mathrm{weighted}} =
\frac{1}{2}\sum_{ij} w_{ij}(X_{ij} - (AS)_{ij})^2
\]

This naturally handles:

- negative flux values
- heteroskedastic noise
- missing data via `w_{ij} = 0`

### 3. Robust loss

For outliers or imperfect sky subtraction, a Student-t-like or Huber-like loss
can be added later.

This is a natural extension point, but not required for a first v2
implementation.

## Regularization

Regularization should remain modular.

### Spectral smoothness

For the component spectra:

\[
\Omega_{\mathrm{smooth}}(S) =
\sum_{r=1}^{k}\sum_{\lambda=2}^{p}(S_{r,\lambda} - S_{r,\lambda-1})^2
\]

This corresponds to the current `lambda_smooth` concept.

### Optional sparsity on abundances

\[
\Omega_{\mathrm{sparse}}(A) = \sum_{i,r} A_{ir}
\]

Useful if later we want simpler or more localized mixing patterns.

### Optional spatial regularization

For IFU data:

\[
\Omega_{\mathrm{spatial}}(A) =
\sum_{r=1}^{k}\sum_{(i,j)\in E}(A_{ir} - A_{jr})^2
\]

where `E` is a graph of neighboring spaxels.

This should be treated as an IFU-specific extension, not as a required part of
the base v2 model.

## Full objective

The general v2 objective becomes:

\[
\mathcal{L}(U_A, U_S) =
\mathcal{L}_{\mathrm{obs}}(X, AS; W)
+ \lambda_{\mathrm{smooth}} \Omega_{\mathrm{smooth}}(S)
+ \lambda_{\mathrm{sparse}} \Omega_{\mathrm{sparse}}(A)
+ \lambda_{\mathrm{spatial}} \Omega_{\mathrm{spatial}}(A)
\]

with

\[
A = \operatorname{softplus}(U_A) + \epsilon
\]

\[
S = \operatorname{softplus}(U_S) + \epsilon
\]

This is effectively a MAP-style estimation framework.

## Relation to Nearly-NMF

Nearly-NMF fits naturally into this broader formulation as one special case:

- `A` and `S` constrained non-negative
- `X` allowed to be signed
- weighted squared-error type objective
- specialized multiplicative updates

That means Nearly-NMF is best understood as:

- a useful reference method
- a possible optional solver
- not the full conceptual foundation of the package

For `SpectralUnmix`, the broader latent-variable framing is the better anchor.

## API direction

A backward-compatible future interface could look like:

```r
spectral_unmix(
  x,
  k,
  obs_model = c("gaussian", "weighted_gaussian", "student_t"),
  weights = NULL,
  positivity = c("softplus", "clamp"),
  lambda_smooth = 0.01,
  lambda_sparse = 0,
  lambda_spatial = 0,
  solver = c("adam"),
  nstart = 1,
  seed = NULL
)
```

Recommended defaults:

- `obs_model = "gaussian"` when `weights` is `NULL`
- `obs_model = "weighted_gaussian"` when `weights` are supplied
- `positivity = "softplus"` for the new implementation
- keep `clamp` as a legacy or compatibility option initially

## Prediction for new data

Prediction should be framed as inference on new abundances with fixed learned
components:

\[
X_{\mathrm{new}} \approx A_{\mathrm{new}} S
\]

That means:

- keep `S` fixed
- optimize only `A_new`
- use the same positivity parameterization
- use the same observation model and optional weights

This matches the current package structure well.

## Minimal implementation plan

### Phase 1

- keep current public API stable
- add `obs_model`
- add optional `weights`
- add `positivity = "softplus"`
- implement weighted Gaussian loss in the torch optimizer

### Phase 2

- make `softplus` the preferred positivity mechanism
- retain clamping only for backward compatibility
- extend `predict()` to weighted fitting for new data

### Phase 3

- optional robust observation model
- optional IFU-specific spatial regularization
- optional reference solver such as Nearly-NMF or other NMF-style updates for
  benchmark comparisons

## What should stay out of scope for now

To keep the package coherent, the following should not be mixed into the first
v2 implementation:

- a full Bayesian sampler
- highly instrument-specific likelihoods
- multiple solver families exposed at once without clear motivation
- public API churn that breaks the simple current workflow

## Bottom line

The most useful generalization for `SpectralUnmix` is:

- positive latent factors
- signed noisy observations
- explicit observation model
- modular regularization
- torch autodiff optimization

This gives a single framework that covers:

- classical NMF-like behavior
- negative-valued data
- inverse-variance weighting
- smooth spectral components
- future IFU-aware penalties

and does so without forcing the package to revolve around one specific update
rule from the literature.
