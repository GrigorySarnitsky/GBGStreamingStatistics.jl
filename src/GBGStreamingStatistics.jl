# GBGStreamingStatistics — umbrella module: it sets the include order and assembles the public API.
#
# The conceptual model (empirical measure → mergeable summary → post-hoc-smoothed estimator; samples
# in, statistics out; estimand-agnostic accumulation, post-hoc smoothing), the invariants, and the
# estimand-first writing/docstring conventions are not duplicated here — they live in
# `docs/ARCHITECTURE.md` (the rationale) and `CLAUDE.md` (the operational rules).


"""
    GBGStreamingStatistics

A backend-agnostic library for **online nonparametric statistics**: it estimates marginal densities,
conditional expectations, conditional and two-time transition densities, and standard functionals of
these, from a stream of samples that is never stored. Each sample is observed once, folded into a
fixed-size summary, and discarded.

The summary is accumulated online and is estimand-agnostic; every reported quantity is a post-hoc
smoothing functional of that one summary — so there is no storage of the underlying state and no
recompute-from-raw reruns. The library is domain-neutral and model-neutral: it returns estimators and
the diagnostics that *test* a model (the Pawula ratio, a Markov / Chapman–Kolmogorov test), never the
model itself.

The package provides the mergeable accumulators and histogram summaries (`AbstractAccumulator`,
`observe!`, `observe_pairs!`, `merge!`, `scale!`, `bin_index`, `Histogram1D`, `Histogram2D` — the last
two with a `count = Float64`
cone variant for hyperstatistics / forgetting), the marginal density (`density` and its variability band
`density_standard_error`, with `gaussian_kernel` / `epanechnikov_kernel`, the half-line `gamma_kernel` and bounded `beta_kernel`, and the parametric starts `normal_start` / `lognormal_start` / `gamma_start` / `exponential_start` / `inverse_gaussian_start` — both the asymmetric kernels and the starts extensible through the exported `AsymmetricKernel` / `ParametricStart` constructors), conditional moments
(`conditional_mean`, `conditional_variance`), conditional and transition densities
(`conditional_density`, `transition_density`), and bandwidth selection (`bandwidth`) — plus the
power-sum quadrature readouts: the recovered marginal moments (`mean`, `var`, `std`, `central_moment`),
the expectation of a smooth functional of ``X`` (`expectation`, with its sampling
`expectation_standard_error`), and the pushforward density of a transformed variable ``g(X)``
(`pushforward_density`).
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
  ParametricStart,
  normal_start,
  lognormal_start,
  gamma_start,
  exponential_start,
  inverse_gaussian_start


using .ConditionalExpectations:
  conditional_mean,
  conditional_variance


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
       conditional_density,
       transition_density,
       bandwidth,
       expectation,
       expectation_standard_error,
       pushforward_density


end # module GBGStreamingStatistics
