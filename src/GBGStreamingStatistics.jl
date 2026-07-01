# GBGStreamingStatistics — umbrella module: it sets the include order and assembles the public API.
#
# The conceptual model (empirical measure → mergeable summary → post-hoc-smoothed estimator; samples
# in, statistics out; estimand-agnostic accumulation, post-hoc smoothing), the invariants, and the
# estimand-first writing/docstring conventions are not duplicated here — they live in
# `docs/ARCHITECTURE.md` (the rationale) and `CLAUDE.md` (the operational rules).


"""
    GBGStreamingStatistics

A backend-agnostic library for **online nonparametric statistics**: it estimates marginal densities,
conditional expectations, conditional and two-time transition densities, and functionals of these
(expectations, moments, and pushforward densities), from a stream of samples that is never stored. Each
sample is observed once, folded into a fixed-size summary, and discarded.

The summary is accumulated online and is estimand-agnostic; every reported quantity is a post-hoc
smoothing functional of that one summary — so there is no storage of the underlying state and no
recompute-from-raw reruns. The library is domain-neutral and model-neutral: no stationarity, Markov, or
parametric structure is a premise — such structure is only ever *tested*, never assumed.

The package provides the mergeable **summaries** (`Histogram1D` / `Histogram2D`, with `observe!` /
`observe_pairs!` / `merge!` / `scale!` — a `count = Float64` cone variant for hyperstatistics and
forgetting); the marginal **`density`** (with its variability band, symmetric or bounded-support kernels,
an adaptive bandwidth, and a parametric start); **conditional expectations** (`conditional_mean`,
`conditional_variance`, and the general `conditional_expectation`) and the conditional / transition
**densities** (`conditional_density`, `transition_density`); **`bandwidth`** selection; and the power-sum
quadrature **readouts** — the moments (`mean` / `var` / `std` / `central_moment`), `expectation`, and
`pushforward_density`. Every estimator has a matching `*_standard_error`, with an estimate-plus-error
`*_with_standard_error` pairing; the API reference documents the full surface.
"""
module GBGStreamingStatistics


include("Accumulators.jl")
include("Histograms.jl")
include("Kernels.jl")
include("Expectations.jl")
include("LocalPolynomial.jl")
include("MarginalDensity.jl")
include("ConditionalExpectations.jl")
include("ConditionalDensity.jl")
include("Bandwidth.jl")
# include("IncrementMoments.jl")         # [archived] the §4.4 SDE-reading functionals — see archive/
# include("CharacteristicFunctions.jl")  # [deferred] the empirical-characteristic-function summary
# include("ProjectionInference.jl")      # [deferred] the SFI / ULI projection road


using .Accumulators:
  AbstractAccumulator,
  observe!,
  merge!,
  scale!


using .Histograms:
  Histogram1D,
  Histogram2D,
  bin_index,
  observe_pairs!


using .Kernels:
  gaussian_kernel,
  epanechnikov_kernel,
  AsymmetricKernel,
  gamma_kernel,
  beta_kernel


using .MarginalDensity:
  density,
  density_standard_error,
  density_with_standard_error,
  ParametricStart,
  normal_start,
  lognormal_start,
  gamma_start,
  exponential_start,
  inverse_gaussian_start


using .ConditionalExpectations:
  conditional_mean,
  conditional_variance,
  conditional_expectation,
  conditional_expectation_standard_error,
  conditional_expectation_with_standard_error,
  conditional_pushforward_density


using .ConditionalDensity:
  conditional_density,
  transition_density


using .Bandwidth:
  bandwidth


using .Expectations:
  mean,
  var,
  std,
  central_moment,
  expectation,
  expectation_standard_error,
  expectation_with_standard_error,
  pushforward_density


export AbstractAccumulator,
       Histogram1D,
       Histogram2D,
       observe!,
       observe_pairs!,
       merge!,
       scale!,
       bin_index,
       mean,
       var,
       std,
       central_moment,
       density,
       density_standard_error,
       density_with_standard_error,
       gaussian_kernel,
       epanechnikov_kernel,
       AsymmetricKernel,
       gamma_kernel,
       beta_kernel,
       ParametricStart,
       normal_start,
       lognormal_start,
       gamma_start,
       exponential_start,
       inverse_gaussian_start,
       conditional_mean,
       conditional_variance,
       conditional_expectation,
       conditional_expectation_standard_error,
       conditional_expectation_with_standard_error,
       conditional_pushforward_density,
       conditional_density,
       transition_density,
       bandwidth,
       expectation,
       expectation_standard_error,
       expectation_with_standard_error,
       pushforward_density


end # module GBGStreamingStatistics
