# GBGStreamingStatistics.jl

`GBGStreamingStatistics.jl` is a backend-agnostic library for online nonparametric statistics. It
estimates marginal densities, conditional expectations, and conditional and two-time transition
densities, from a stream of samples that is never stored — each sample is observed once,
folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin power-sums about the bin center) is accumulated online and is *estimand-agnostic*:
every reported quantity is a post-hoc smoothing functional of that one summary, so the samples themselves are
never stored and never re-read. It is a mergeable monoid — in fact a
cone (combine across threads, devices, or time-chunks by addition, or weight and superpose for
hyperstatistics — a post-hoc mixture over summaries) — and targets both CPU and GPU array backends on a single shared-memory node. The library is domain-neutral and model-neutral — it assumes no
statistical model (no stationarity, Markov, diffusion, or parametric premise): it consumes a stream
of values (optionally lagged pairs) and returns estimators.

**When it fits.** It is built for settings where the data is a stream and storing or re-reading it is
the bottleneck:

- **inside a time-marching solver** — fold each step's samples into the summary and discard, with no
  trajectory storage and no second pass;
- **one pass, many views** — every bandwidth, estimand, and conditioning value is computed post hoc
  from the one summary, so exploring them never re-reads the raw stream;
- **parallel and chunked** — mergeable summaries combine across threads, GPUs, and time windows;
- **memory-tight or GPU-resident** — the summary is fixed-size and the accumulation hot path is
  allocation-free and runs on the GPU.

It suits, for example, turbulence (conditional and Lagrangian statistics, transition densities)
and molecular dynamics (the van Hove and intermediate scattering functions are
two-time conditional densities) — but nothing in the library is specific to a domain.

!!! note "Status"
    The summaries [`Histogram1D`](@ref) and [`Histogram2D`](@ref) ([`observe!`](@ref) /
    [`observe_pairs!`](@ref) / `merge!` / [`scale!`](@ref), with a floating-point-count variant that scales post hoc for
    hyperstatistics or the exponential forgetting of older samples), the marginal [`density`](@ref), global [`bandwidth`](@ref) selection,
    the conditional mean and variance ([`conditional_mean`](@ref), [`conditional_variance`](@ref)),
    and the conditional and transition densities ([`conditional_density`](@ref), [`transition_density`](@ref))
    are implemented and tested on CPU and GPU (analytic oracles, mergeability, allocation, GPU results
    matching the CPU's), including the asymmetric kernels for bounded or positive data, the local
    bandwidth ``h^\star(x)`` fitted to the response at each query point, and the Fan–Yao conditional
    variance — plus the estimators read by quadrature over the stored power-sums ([`expectation`](@ref),
    its [`expectation_standard_error`](@ref), and [`pushforward_density`](@ref)). Drift and diffusion
    (the first two Kramers–Moyal coefficients) were built and tested, then archived, restorable when a
    consumer (a package built on this one) needs them; the higher coefficients and the Pawula / Markov
    diagnostics are forthcoming. See the [Architecture](architecture.md) for the design and roadmap.

## Installation

Install it by name from the GBG packages' own registry,
[GBGRegistry](https://github.com/GrigorySarnitsky/GBGRegistry), added after `General`:

```julia
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(url = "https://github.com/GrigorySarnitsky/GBGRegistry")
Pkg.add("GBGStreamingStatistics")
```

or at the package prompt:

```
pkg> registry add General https://github.com/GrigorySarnitsky/GBGRegistry
pkg> add GBGStreamingStatistics
```

Then you can use it as usual:

```julia
using GBGStreamingStatistics
```

New here? Start with the [Getting started](tutorial.md) tour — the intended *samples in,
statistics out* workflow.

## Main components

The package splits into a *universal lower half* (the summaries and how they accumulate, shared
verbatim by every consumer) and an *estimator upper half* (the estimands):

- [`Accumulators`](api/Accumulators.md): the foundation — [`AbstractAccumulator`](@ref), the interface of every
  mergeable online summary of the [empirical measure](https://en.wikipedia.org/wiki/Empirical_measure)
  ``\mathbb{P}_n``, which gives a set ``A`` the share ``\mathbb{P}_n(A) = \#\{i : X_i \in A\}/n`` of the samples lying in it (you `observe!` samples into it and `merge!` two together).
- [`Histograms`](api/Histograms.md): the summaries — `Histogram1D` and `Histogram2D`, both
  `<: AbstractAccumulator`. Each is a finely binned empirical measure holding per-bin counts and
  power-sums about the bin center; accumulation runs on CPU or GPU. It also offers `observe_pairs!`
  for two-time statistics and, when `HDF5` is loaded (a package extension),
  [`write_histogram`](@ref) / [`read_histogram`](@ref), persisting a summary to an HDF5 group.
- [`Kernels`](api/Kernels.md): the shared smoothing kernels, which `density`, the regression fits, the
  conditional density and the bandwidth selectors use. They are symmetric
  ([`gaussian_kernel`](@ref), [`epanechnikov_kernel`](@ref), or any callable kernel with ``\int K = 1``)
  or asymmetric for bounded/positive support ([`gamma_kernel`](@ref), [`beta_kernel`](@ref), or a custom
  [`AsymmetricKernel`](@ref)).
- [`Expectations`](api/Expectations.md): the power-sum quadrature readouts (expectations computed by
  contracting the stored per-bin power-sums with a polynomial fit of the functional) — the expectation of a smooth
  functional `expectation` (``\mathbb{E}[g(X)]``), its sampling `expectation_standard_error` (paired
  by [`expectation_with_standard_error`](@ref)), and the pushforward density
  `pushforward_density` (``f_{g(X)}``) of a transformed variable. The moments (`mean` / `var` / `std` /
  [`central_moment`](@ref)) are the monomial case ``g = x^p``, one per order.
- [`MarginalDensity`](api/MarginalDensity.md): the marginal density `density` (``f_X``) offers, in one
  call, the binned KDE, the local-log-polynomial estimator (`degree`), the Abramson variable
  bandwidth (`adaptive`, for tails), the asymmetric `gamma_kernel` / `beta_kernel` for bounded/positive
  support, the parametric `start` (a moment-fitted model corrected by a KDE; Hjort & Glad 1995), and exact
  `normalize` — plus [`density_standard_error`](@ref), the pointwise variability band, and
  [`density_with_standard_error`](@ref), the estimate and its band together.
- [`ConditionalExpectations`](api/ConditionalExpectations.md): the conditional expectation
  [`conditional_expectation`](@ref) (``\mathbb{E}[g(Y)\mid X = x]``) of a smooth response
  functional. The conditional mean `conditional_mean` (``m(x)``) and variance
  `conditional_variance` (``\sigma^2(x)``) sit beside it — the mean as its ``g = \mathrm{identity}`` case,
  the variance's plug-in method built from its ``y^2`` case. It also offers
  [`conditional_expectation_standard_error`](@ref) (paired by
  [`conditional_expectation_with_standard_error`](@ref)) and the conditional pushforward density
  [`conditional_pushforward_density`](@ref) (``f_{g(Y)\mid X}``), each with its standard error. All are
  local-polynomial regressions in ``x`` of a ``y``-quadrature response (per-bin expectations from
  the stored response power-sums).
- [`ConditionalDensity`](api/ConditionalDensity.md): the conditional density `conditional_density`
  (``f_{Y\mid X}``) and two-time transition density `transition_density`
  (``p_\tau(y\mid x)``, the density of the value a lag later) are double-kernel estimators with
  standard errors and paired `_with_standard_error` functions.
- [`Bandwidth`](api/Bandwidth.md): `bandwidth` applies the default kernel- and degree-aware
  normal-reference rule (Silverman's is the Gaussian, degree-0 case), a data-driven plug-in
  (e.g. Sheather–Jones, `rule = :sj`), or, for paired samples, the local plug-in
  ``h^\star(x)`` (`rule = :local`).

## Documentation structure

### Start here

- [Getting started](tutorial.md): the guided tour of the *samples in, statistics out* workflow, from a
  summary to a smoothed estimate and its standard error.
- [Notation](notation.md): the single reference for every symbol the docstrings use.

### API reference

```@contents
Pages = [
  "api.md",
  "api/GBGStreamingStatistics.md",
  "api/Accumulators.md",
  "api/Histograms.md",
  "api/Kernels.md",
  "api/Expectations.md",
  "api/MarginalDensity.md",
  "api/ConditionalExpectations.md",
  "api/ConditionalDensity.md",
  "api/Bandwidth.md",
]
```

### Design & internals

- [Architecture](architecture.md): why the package is built as it is — the empirical-measure spine,
  estimand-agnostic accumulation with post-hoc smoothing, the invariants, and the roadmap.
- [Formalization](formalization.md): the typed statements behind the prose — each estimator's
  estimand, estimator, and consistency claim, linked to the tests.
- [Testing](testing.md): the analytic oracles (Gaussian, Ornstein–Uhlenbeck) and the per-operation
  coverage contract.
- [Performance](performance.md): the performance model — the accumulation hot loop as the single
  critical path.
