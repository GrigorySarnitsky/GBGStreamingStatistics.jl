# GBGStreamingStatistics.jl

`GBGStreamingStatistics.jl` is a backend-agnostic library for **online nonparametric statistics**. It
estimates marginal densities, conditional expectations, and conditional and two-time transition
densities, from a **stream of samples that is never stored** — each sample is observed once,
folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin moments) is accumulated online and is *estimand-agnostic*:
every reported quantity is a post-hoc smoothing functional of that one summary, so there is no storage
of the underlying state and no recompute-from-raw reruns. It is a **mergeable monoid** — in fact a
**cone** (combine across threads, devices, or time-chunks by addition, or weight and superpose for
hyperstatistics) — and targets both **CPU and GPU** array backends on a single shared-memory node. The library is **domain-neutral and model-neutral** — it assumes no
statistical model (no stationarity, Markov, diffusion, or parametric premise): it consumes a stream
of values (optionally lagged pairs) and returns estimators.

**When it fits.** It is built for settings where the data is a stream and storing or re-reading it is
the bottleneck:

- **inside a time-marching solver** — fold each step's samples into the summary and discard, with no
  trajectory storage and no second pass;
- **one pass, many views** — every bandwidth, channel, and conditioning value is computed post hoc
  from the one summary, so exploring them never re-reads the raw stream;
- **parallel and chunked** — mergeable summaries combine across threads, GPUs, and time windows;
- **memory-tight or GPU-resident** — the summary is fixed-size and the accumulation hot path is
  allocation-free and runs on the GPU.

It is used, for example, in turbulence (conditional and Lagrangian statistics, transition densities,
drift and diffusion) and molecular dynamics (the van Hove and intermediate scattering functions are
two-time conditional densities) — but nothing in the library is specific to a domain.

!!! note "Status"
    The summaries `Histogram1D` and `Histogram2D` (power-sums about the bin center, `observe!` /
    `observe_pairs!` / `merge!` / `scale!` / `bin_index`, with a Float-count cone variant for
    hyperstatistics / forgetting), the marginal [`density`](@ref), global [`bandwidth`](@ref) selection,
    and the conditional mean and variance ([`conditional_mean`](@ref), [`conditional_variance`](@ref)),
    the conditional and transition densities ([`conditional_density`](@ref), [`transition_density`](@ref))
    are implemented and tested on **CPU and GPU** (Gaussian, bivariate-Gaussian, and Ornstein–Uhlenbeck
    oracles, mergeability, allocation, AMD-ROCm parity), including the asymmetric bounded/positive-support
    kernels, the response-aware local bandwidth ``h^\star(x)``, and the Fan–Yao conditional variance —
    plus the power-sum quadrature readouts ([`expectation`](@ref) of a smooth functional, its
    [`expectation_standard_error`](@ref), and the [`pushforward_density`](@ref) of a transformed
    variable). The Kramers–Moyal / SDE-reading functionals (drift and diffusion), the higher
    Kramers–Moyal coefficients, and the Pawula / Markov diagnostics are not currently provided but are
    planned; see the [Architecture](architecture.md) for the design and roadmap.

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

New here? Start with the **[Getting started](tutorial.md)** tour — the intended *samples in,
statistics out* workflow.

## Main components

The package splits into a **universal lower half** (the summaries and how they accumulate, shared
verbatim by every consumer) and an **estimator upper half** (the estimands):

- [`Accumulators`](api/Accumulators.md): the foundation — `AbstractAccumulator`, the interface of every
  mergeable online **summary** of the [empirical measure](https://en.wikipedia.org/wiki/Empirical_measure)
  ``\hat\mu_n = \frac1n\sum_i \delta_{X_i}`` (you `observe!` samples into it and `merge!` two together).
  The concrete summaries live in `Histograms`.
- [`Histograms`](api/Histograms.md): the summaries — `Histogram1D` and `Histogram2D`, both
  `<: AbstractAccumulator`. Each is a finely **binned empirical measure** — the *fine-grained measure*
  that post-hoc smoothing later coarse-grains — holding per-bin counts and power-sums about the bin
  center; accumulation runs on **CPU or GPU**, and summaries combine with `merge!`. With `bin_index`
  and `observe_pairs!` for two-time statistics.
- [`Kernels`](api/Kernels.md): the shared smoothing kernels — symmetric (`gaussian_kernel`, or any
  callable with ``\int K = 1``) and asymmetric bounded/positive-support kernels (`gamma_kernel`, via the
  `AsymmetricKernel` constructor) — used by `density`, the regression fits, the double-kernel conditional
  density, and the bandwidth rules.
- [`Expectations`](api/Expectations.md): the power-sum quadrature readouts — the expectation of a smooth
  functional `expectation` ``\equiv \mathbb{E}[g(X)]``, its sampling `expectation_standard_error` (paired
  with the estimate by `expectation_with_standard_error`), and the pushforward density
  `pushforward_density` ``\equiv f_{g(X)}`` of a transformed variable. Read off the stored power-sums; the
  recovered moments (`mean` / `var` / `std` / `central_moment`) are the ``g = x^p`` case.
- [`MarginalDensity`](api/MarginalDensity.md): the marginal density `density` ``\equiv f_X`` — one
  call over the binned KDE, the local-log-polynomial (`degree`), and the Abramson variable bandwidth
  (`adaptive`, for tails), the asymmetric `gamma_kernel` / `beta_kernel` for bounded/positive support,
  the parametric `start` (a moment-fitted model corrected by a KDE; Hjort & Glad 1995), and exact
  `normalize` — plus `density_standard_error`, the pointwise variability band (tracking those same
  options), and `density_with_standard_error`, the estimate and its band from one option specification.
- [`ConditionalExpectations`](api/ConditionalExpectations.md): the conditional expectation
  `conditional_expectation` ``\equiv \mathbb{E}[g(Y)\mid x]`` of a response functional, with the conditional
  mean `conditional_mean` ``\equiv m(x)`` and variance `conditional_variance` ``\equiv \sigma^2(x)`` as its
  ``g = \mathrm{identity}`` / ``y^2`` cases, the sampling `conditional_expectation_standard_error` (paired
  with the estimate by `conditional_expectation_with_standard_error`), and the
  conditional pushforward density `conditional_pushforward_density` ``\equiv f_{g(Y)\mid X}``
  (local-polynomial regression in ``x`` of a ``y``-quadrature response).
- [`ConditionalDensity`](api/ConditionalDensity.md): the conditional density `conditional_density`
  ``\equiv f_{Y\mid X}`` and two-time transition density `transition_density`
  ``\equiv p_\tau(y\mid x)`` (double-kernel).
- [`Bandwidth`](api/Bandwidth.md): bandwidth selection — the default kernel- and degree-aware
  normal-reference rule (Silverman the Gaussian/degree-0 case), plus data-driven plug-ins
  (e.g. Sheather–Jones, `rule = :sj`) and the joint-summary response-aware local plug-in ``h^\star(x)``
  (`rule = :local`), among others.

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
