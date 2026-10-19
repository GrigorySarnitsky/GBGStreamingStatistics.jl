# GBGStreamingStatistics — umbrella module: it sets the include order and assembles the public API.
#
# The conceptual model (empirical measure → mergeable summary → post-hoc-smoothed estimator; samples
# in, statistics out; estimand-agnostic accumulation, post-hoc smoothing), the invariants, and the
# estimand-first writing/docstring conventions are not duplicated here — they live in
# `docs/ARCHITECTURE.md` (the rationale) and `CLAUDE.md` (the operational rules).


"""
    GBGStreamingStatistics

A library for online nonparametric statistics, estimating densities, conditional
expectations and functionals of these from a stream of samples: each sample is observed once, folded
into a fixed-size summary, and discarded. It runs on CPUs and GPUs.

The summary records per bin the count and the power sums of its samples about the bin center. It is
estimand-agnostic: every reported quantity is a post-hoc smoothing functional of it, and no stationarity,
Markov, or parametric structure is a premise.

It offers these operation families:

- Summaries — [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) /
  [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D), filled by
  [`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!) /
  [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!), added by
  [`merge!`](@ref Base.merge!) and, with floating-point counts (`count_type = Float64`), reweighted by
  [`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!) for mixtures weighted after accumulation
  (hyperstatistics) and exponential forgetting.
- [Smoothing kernels](@ref "Kernels").
- Marginal density — [`density`](@ref GBGStreamingStatistics.MarginalDensity.density), optionally from a
  [parametric start](@ref "Parametric starts").
- Conditional expectations —
  [`conditional_expectation`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_expectation),
  [`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean),
  [`conditional_variance`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_variance),
  [`conditional_pushforward_density`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_pushforward_density).
- Conditional and transition densities —
  [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density) /
  [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density).
- Expectations by quadrature over the power sums —
  [`expectation`](@ref GBGStreamingStatistics.Expectations.expectation),
  [`central_moment`](@ref GBGStreamingStatistics.Expectations.central_moment),
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
