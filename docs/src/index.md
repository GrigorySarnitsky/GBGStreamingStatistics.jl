# GBGStreamingStatistics.jl

`GBGStreamingStatistics.jl` is a backend-agnostic library for online nonparametric statistics. It
estimates marginal densities, conditional expectations, and conditional and two-time transition
densities, from a stream of samples that is never stored — each sample is observed once,
folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin power-sums about the bin center) is accumulated online and is *estimand-agnostic*:
every reported quantity is a post-hoc smoothing functional of that one summary, so the underlying state is
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

It is used, for example, in turbulence (conditional and Lagrangian statistics, transition densities,
drift and diffusion) and molecular dynamics (the van Hove and intermediate scattering functions are
two-time conditional densities) — but nothing in the library is specific to a domain.

!!! note "Status"
    The summaries [`Histogram1D`](@ref) and [`Histogram2D`](@ref) (power-sums about the bin center, [`observe!`](@ref) /
    [`observe_pairs!`](@ref) / `merge!` / [`scale!`](@ref), with a Float-count cone variant, re-weightable post hoc for
    hyperstatistics / forgetting), the marginal [`density`](@ref), global [`bandwidth`](@ref) selection,
    the conditional mean and variance ([`conditional_mean`](@ref), [`conditional_variance`](@ref)),
    and the conditional and transition densities ([`conditional_density`](@ref), [`transition_density`](@ref))
    are implemented and tested on CPU and GPU (Gaussian, bivariate-Gaussian, and Ornstein–Uhlenbeck
    oracles, mergeability, allocation, AMD-ROCm parity), including the asymmetric bounded-support
    kernels, the response-aware local bandwidth ``h^\star(x)``, and the Fan–Yao conditional variance —
    plus the readouts by quadrature over the stored power-sums ([`expectation`](@ref) of a smooth
    functional, its [`expectation_standard_error`](@ref), and the [`pushforward_density`](@ref) of a
    transformed variable). The Kramers–Moyal functionals (drift, diffusion, the higher coefficients) and
    the Pawula / Markov diagnostics are archived — built and tested, then moved out to keep the core
    focused, and restorable when a consumer (a package built on this one) needs them; see the [Architecture](architecture.md) for the design and roadmap.

## Installation

This is an unregistered package, so it is not listed in the `General` Julia registry. To install it, run:

```julia
]add https://github.com/GrigorySarnitsky/GBGStreamingStatistics.jl
```

or

```julia
using Pkg

Pkg.add(url="https://github.com/GrigorySarnitsky/GBGStreamingStatistics.jl")
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
  ``\mathbb{P}_n``, which gives a set ``A`` the share ``\mathbb{P}_n(A) = \#\{i : X_i \in A\}/n`` of the samples (you `observe!` samples into it and `merge!` two together).
- [`Histograms`](api/Histograms.md): the summaries — `Histogram1D` and `Histogram2D`, both
  `<: AbstractAccumulator`. Each is a finely binned empirical measure holding per-bin counts and
  power-sums about the bin center; accumulation runs on CPU or GPU, and summaries combine with
  `merge!`. It also offers `observe_pairs!` for two-time statistics and, when `HDF5` is loaded (a
  package extension), [`write_histogram`](@ref) / [`read_histogram`](@ref), which persist a summary to an HDF5 group.
- [`Kernels`](api/Kernels.md): the shared smoothing kernels — symmetric ([`gaussian_kernel`](@ref),
  [`epanechnikov_kernel`](@ref), or any callable with ``\int K = 1``) and asymmetric bounded/positive-support
  kernels ([`gamma_kernel`](@ref), [`beta_kernel`](@ref), or a custom [`AsymmetricKernel`](@ref)) — used by `density`, the
  regression fits, the double-kernel conditional
  density, and the bandwidth rules.
- [`Expectations`](api/Expectations.md): the power-sum quadrature readouts (expectations evaluated by
  contracting the stored per-bin power-sums with a polynomial fit of the functional) — the expectation of a smooth
  functional `expectation` (``\mathbb{E}[g(X)]``), its sampling `expectation_standard_error` (paired
  with the estimate by [`expectation_with_standard_error`](@ref)), and the pushforward density
  `pushforward_density` (``f_{g(X)}``) of a transformed variable. All are read off the stored power-sums;
  the recovered moments (`mean` / `var` / `std` / [`central_moment`](@ref)) are the ``g = x^p`` case.
- [`MarginalDensity`](api/MarginalDensity.md): the marginal density `density` (``f_X``) — one
  call over the binned KDE, the local-log-polynomial estimator (`degree`), the Abramson variable bandwidth
  (`adaptive`, for tails), the asymmetric `gamma_kernel` / `beta_kernel` for bounded/positive support,
  the parametric `start` (a moment-fitted model corrected by a KDE; Hjort & Glad 1995), and exact
  `normalize` — plus [`density_standard_error`](@ref), the pointwise variability band, and
  [`density_with_standard_error`](@ref), the estimate and its band together.
- [`ConditionalExpectations`](api/ConditionalExpectations.md): the conditional expectation
  [`conditional_expectation`](@ref) (``\mathbb{E}[g(Y)\mid X = x]``) of the response, with the conditional mean `conditional_mean` (``m(x)``) and variance
  `conditional_variance` (``\sigma^2(x)``, with its residual-based standard error) as its
  ``g = \mathrm{identity}`` / ``y^2`` cases, the sampling [`conditional_expectation_standard_error`](@ref) (paired
  with the estimate by [`conditional_expectation_with_standard_error`](@ref)), and the conditional pushforward
  density [`conditional_pushforward_density`](@ref) (``f_{g(Y)\mid X}``) — each with its standard error. All are
  local-polynomial regressions in ``x`` of a ``y``-quadrature response (the per-bin expectation of the
  functional, read off the stored response power-sums).
- [`ConditionalDensity`](api/ConditionalDensity.md): the conditional density `conditional_density`
  (``f_{Y\mid X}``) and two-time transition density `transition_density`
  (``p_\tau(y\mid x)``; double-kernel), each with its pointwise standard error and a
  `_with_standard_error` reader returning both together.
- [`Bandwidth`](api/Bandwidth.md): bandwidth selection — the default kernel- and degree-aware
  normal-reference rule (Silverman's rule is the Gaussian, degree-0 case), plus data-driven plug-ins
  (e.g. Sheather–Jones, `rule = :sj`) and the joint-summary response-aware local plug-in ``h^\star(x)``
  (`rule = :local`).

## Documentation structure

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
