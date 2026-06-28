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

The package provides the mergeable accumulators and histogram summaries (`observe!`, `observe_pairs!`,
`merge!`, `bin_index`, `Histogram1D`, `Histogram2D`) with their recovered marginal moments (`mean`,
`var`, `std`, `central_moment`), the marginal density (`density` and its variability band
`density_standard_error`, with `gaussian_kernel` / `epanechnikov_kernel`, the bounded-support `gamma_kernel` / `beta_kernel`, and the parametric starts `normal_start` / `lognormal_start` / `gamma_start` / `exponential_start` / `inverse_gaussian_start`), conditional moments
(`conditional_mean`, `conditional_variance`), conditional and transition densities
(`conditional_density`, `transition_density`), the Kramers–Moyal functionals (`increment_moment`,
`kramers_moyal`, `drift`, `diffusion`, with the `pawula_quotient` and `markov_test` diagnostics
forthcoming), and bandwidth selection (`bandwidth`).
"""
module GBGStreamingStatistics


include("Accumulators.jl")
include("Histograms.jl")
include("DensityEstimation.jl")
include("Regression.jl")
include("ConditionalDensity.jl")
include("MarkovInference.jl")
include("Bandwidth.jl")
# include("CharacteristicFunctions.jl")  # [deferred] the empirical-characteristic-function summary
# include("ProjectionInference.jl")      # [deferred] the SFI / ULI projection road


using .Accumulators:
  AbstractAccumulator,
  observe!,
  merge!


using .Histograms:
  Histogram1D,
  Histogram2D,
  bin_index,
  observe_pairs!,
  mean,
  var,
  std,
  central_moment


using .DensityEstimation:
  density,
  density_standard_error,
  gaussian_kernel,
  epanechnikov_kernel,
  gamma_kernel,
  beta_kernel,
  normal_start,
  lognormal_start,
  gamma_start,
  exponential_start,
  inverse_gaussian_start


using .Regression:
  conditional_mean,
  conditional_variance


using .ConditionalDensity:
  conditional_density,
  transition_density


using .MarkovInference:
  increment_moment,
  kramers_moyal,
  drift,
  diffusion,
  pawula_quotient,
  markov_test


using .Bandwidth:
  bandwidth


export AbstractAccumulator,
       Histogram1D,
       Histogram2D,
       observe!,
       observe_pairs!,
       merge!,
       bin_index,
       mean,
       var,
       std,
       central_moment,
       density,
       density_standard_error,
       gaussian_kernel,
       epanechnikov_kernel,
       gamma_kernel,
       beta_kernel,
       normal_start,
       lognormal_start,
       gamma_start,
       exponential_start,
       inverse_gaussian_start,
       conditional_mean,
       conditional_variance,
       conditional_density,
       transition_density,
       increment_moment,
       kramers_moyal,
       drift,
       diffusion,
       pawula_quotient,
       markov_test,
       bandwidth


end # module GBGStreamingStatistics
