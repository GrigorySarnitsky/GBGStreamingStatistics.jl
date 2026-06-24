# GBGStreamingStatistics.jl

`GBGStreamingStatistics.jl` is a backend-agnostic library for **online nonparametric statistics**. It
estimates marginal densities, conditional expectations, conditional and two-time transition
densities, and standard functionals of these (Kramers–Moyal coefficients, drift and diffusion, the
Pawula diagnostic), from a **stream of samples that is never stored** — each sample is observed once,
folded into a fixed-size summary, and discarded.

The summary (a fine histogram with per-bin moments) is accumulated online and is *response-agnostic*;
every reported quantity is a post-hoc smoothing functional of that one summary — so there is no
storage of the underlying state and no recompute-from-raw reruns. The library is **domain-neutral and
model-neutral**: it consumes a stream of values (optionally lagged pairs) and returns estimators, with
no concept of velocity, droplet, trajectory, or turbulence, and no assumption that the data is
stationary, Markov, diffusive, or parametric — such structure is only ever *tested*. The
implementation targets both CPU and GPU array backends, on a single shared-memory node.

This package is the statistics sibling of [`GBGPeriodicFields.jl`](https://github.com/GrigorySarnitsky/GBGPeriodicFields.jl): that package is all mathematics and no physics; this one is all statistics and no physics.

!!! note "Status"
    The summaries `Histogram1D` and `Histogram2D` (power-sums about the bin center, `observe!` /
    `observe_pairs!` / `merge!` / `bin_index`), the marginal
    [`density`](@ref GBGStreamingStatistics.DensityEstimation.density) (§4.1), global
    [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) selection, and the conditional mean
    and variance ([`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean),
    `conditional_variance`, §4.2), and the Kramers–Moyal
    [`drift`](@ref GBGStreamingStatistics.MarkovInference.drift) /
    [`diffusion`](@ref GBGStreamingStatistics.MarkovInference.diffusion) (§4.4) are implemented and
    tested on **CPU and GPU** (Gaussian, bivariate-Gaussian, and Ornstein–Uhlenbeck oracles,
    mergeability, allocation, AMD-ROCm parity). The §4.3 conditional and transition densities, the
    response-aware local bandwidth, and the higher Kramers–Moyal / Pawula / Markov diagnostics are
    forthcoming; see the
    [Architecture](architecture.md) for the design and roadmap.

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
  The concrete summaries and their scatter kernels live in `Histograms`.
- [`Histograms`](api/Histograms.md): the summaries — `Histogram1D` and `Histogram2D`, both
  `<: AbstractAccumulator`. Each is a finely **binned empirical measure** — the *fine-grained measure*
  that post-hoc smoothing later coarse-grains — holding per-bin counts and power-sums about the bin
  center, accumulated by one KernelAbstractions + Atomix atomic scatter-add (`merge!` by addition);
  with `bin_index` and `observe_pairs!` for two-time statistics.
- [`DensityEstimation`](api/DensityEstimation.md): the marginal density `density` ``\equiv f_X``
  (binned KDE / WARPing).
- [`Regression`](api/Regression.md): the conditional mean `conditional_mean` ``\equiv m(x)`` and
  variance `conditional_variance` ``\equiv \sigma^2(x)`` (local-linear).
- [`ConditionalDensity`](api/ConditionalDensity.md): the conditional density `conditional_density`
  ``\equiv f_{Y\mid X}`` and two-time transition density `transition_density`
  ``\equiv p_\tau(y\mid x)`` (double-kernel).
- [`MarkovInference`](api/MarkovInference.md): conditional increment moments and their functionals —
  `kramers_moyal` ``\equiv D^{(n)}``, `drift` ``\equiv D^{(1)}``, `diffusion` ``\equiv D^{(2)}``, the
  `pawula_quotient`, and the `markov_test`.
- [`Bandwidth`](api/Bandwidth.md): bandwidth selection — global rules and the response-aware local
  plug-in ``h^\star(x)``.

## Documentation structure

### Design & internals

- [Architecture](architecture.md): why the package is built as it is — the empirical-measure spine,
  response-agnostic accumulation with post-hoc smoothing, the invariants, and the roadmap.
- [Formalization](formalization.md): the typed statements behind the prose — each estimator's
  estimand, estimator, and consistency claim, linked to the tests.
- [Testing](testing.md): the analytic oracles (Gaussian, Ornstein–Uhlenbeck) and the per-operation
  coverage contract.
- [Performance](performance.md): the performance model — the accumulation hot loop as the single
  critical path.

### API reference

```@contents
Pages = [
  "api.md",
  "api/GBGStreamingStatistics.md",
  "api/Accumulators.md",
  "api/Histograms.md",
  "api/DensityEstimation.md",
  "api/Regression.md",
  "api/ConditionalDensity.md",
  "api/MarkovInference.md",
  "api/Bandwidth.md",
]
```
