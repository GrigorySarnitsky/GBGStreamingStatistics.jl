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

On top of the mergeable **summaries** (`Histogram1D` / `Histogram2D`, with a `count = Float64` cone
variant for hyperstatistics and forgetting), it provides — post hoc, from that one summary — marginal
densities, conditional expectations, conditional and two-time transition densities, bandwidth selection,
the power-sum quadrature readouts (moments, expectations, and pushforward densities), and standard-error
bands for the point estimators. See the API reference for the full, module-by-module surface.
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
  conditional_variance_standard_error,
  conditional_variance_with_standard_error,
  conditional_expectation,
  conditional_expectation_standard_error,
  conditional_expectation_with_standard_error,
  conditional_pushforward_density,
  conditional_pushforward_density_standard_error,
  conditional_pushforward_density_with_standard_error


using .ConditionalDensity:
  conditional_density,
  conditional_density_standard_error,
  conditional_density_with_standard_error,
  transition_density,
  transition_density_standard_error,
  transition_density_with_standard_error


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
  pushforward_density,
  pushforward_density_standard_error,
  pushforward_density_with_standard_error


export AbstractAccumulator,
       Histogram1D,
       Histogram2D,
       observe!,
       observe_pairs!,
       merge!,
       scale!,
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
       conditional_variance_standard_error,
       conditional_variance_with_standard_error,
       conditional_expectation,
       conditional_expectation_standard_error,
       conditional_expectation_with_standard_error,
       conditional_pushforward_density,
       conditional_pushforward_density_standard_error,
       conditional_pushforward_density_with_standard_error,
       conditional_density,
       conditional_density_standard_error,
       conditional_density_with_standard_error,
       transition_density,
       transition_density_standard_error,
       transition_density_with_standard_error,
       bandwidth,
       expectation,
       expectation_standard_error,
       expectation_with_standard_error,
       pushforward_density,
       pushforward_density_standard_error,
       pushforward_density_with_standard_error


end # module GBGStreamingStatistics
