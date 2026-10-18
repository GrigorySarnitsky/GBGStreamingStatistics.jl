# GBGStreamingStatistics — umbrella module: it sets the include order and assembles the public API.
#
# The conceptual model (empirical measure → mergeable summary → post-hoc-smoothed estimator; samples
# in, statistics out; estimand-agnostic accumulation, post-hoc smoothing), the invariants, and the
# estimand-first writing/docstring conventions are not duplicated here — they live in
# `docs/ARCHITECTURE.md` (the rationale) and `CLAUDE.md` (the operational rules).


"""
    GBGStreamingStatistics

A backend-agnostic library for online nonparametric statistics, estimating densities, conditional
expectations and functionals of these from a stream of samples that is never stored: each sample is
observed once, folded into a fixed-size summary, and discarded.

The summary is accumulated online and is estimand-agnostic: every reported quantity is a post-hoc
smoothing functional of it, so the underlying state is never stored or re-read; no stationarity, Markov,
or parametric structure is a premise.

Its operation families, each read post hoc from the one summary:

- Summaries and accumulation — [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) /
  [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D): mergeable (two summaries add
  into one), the `count_type = Float64` variant also rescalable for post-hoc mixtures
  (hyperstatistics) and forgetting.
- Marginal density — [`density`](@ref GBGStreamingStatistics.MarginalDensity.density).
- Conditional expectations —
  [`conditional_expectation`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_expectation),
  the conditional mean and variance.
- Conditional and transition densities —
  [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density) /
  [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density).
- Quadrature readouts over the stored power sums —
  [`expectation`](@ref GBGStreamingStatistics.Expectations.expectation), the moments,
  [`pushforward_density`](@ref GBGStreamingStatistics.Expectations.pushforward_density).
- Bandwidth selection — [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth).
- HDF5 persistence — [`write_histogram`](@ref GBGStreamingStatistics.Histograms.write_histogram) /
  [`read_histogram`](@ref GBGStreamingStatistics.Histograms.read_histogram) (needs `HDF5`).
"""
module GBGStreamingStatistics


# `bandwidth` is declared here, ahead of the modules: `density`'s default bandwidth names it
# (`h = bandwidth(summary; kernel, degree)`), while its methods live in `Bandwidth`, included last.
function bandwidth end

include("Accumulators.jl")
include("Histograms.jl")
include("Kernels.jl")
include("Expectations.jl")
include("LocalPolynomial.jl")
include("MarginalDensity.jl")
include("ConditionalExpectations.jl")
include("ConditionalDensity.jl")
include("Bandwidth.jl")
# include("IncrementMoments.jl")         # [archived] the §4.4 SDE-reading functionals — git commit 37006b5, archive/
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
  observe_pairs!,
  write_histogram,
  read_histogram


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
       pushforward_density_with_standard_error,
       write_histogram,
       read_histogram


end # module GBGStreamingStatistics
