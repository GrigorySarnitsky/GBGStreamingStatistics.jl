# Oracle provenance — where each public operation's ground truth comes from.
#
# The suite's bar (CLAUDE.md, Definition of done) is an INDEPENDENT analytic oracle for every estimator:
# a Gaussian or Ornstein–Uhlenbeck closed form, a Monte-Carlo sampling spread, a direct computation on
# the raw data — never only agreement with another operation of this package. The registry is DECLARED to
# GBGCodingTools (RULE.TEST.ORACLE_INDEPENDENT): the engine is the tool's `TestHarness.oracle_registry`.
# It asserts that every public value-returning operation is classified here, that no entry names a
# vanished operation, that every class is one of the five below, that no operation rests on `relational`
# coverage alone, and that every non-structural citation names a test file that exists and mentions the
# operation. A new operation fails this file until someone classifies it — which is the moment to notice
# that it has no independent oracle.
#
#   analytic     a closed form derived independently of the package (the Gaussian pdf, E[cos X] = e^{-1/2},
#                the OU transition law, Silverman's constant)
#   independent  a direct computation written in the test (the moments from the raw samples, the
#                Monte-Carlo spread of an estimate across independent datasets)
#   property     an invariant that pins the operation without computing it
#   relational   compared ONLY against other package operations — legitimate as a second layer, never as
#                the only one
#   structural   returns metadata or storage, not a value
#
# The table stays here, in Julia, so each entry can carry its explanation; the engine is handed the same
# entries in the shape it reads.

using Test
using GBGStreamingStatistics
import GBGCodingTools
include("public_surface.jl")

# operation => (provenance, where the strongest oracle lives)
const ORACLE_PROVENANCE = Dict{Symbol, Tuple{Symbol, String}}(
  # ---- the recovered moments: exact against a direct computation on the raw samples ----
  :mean            => (:independent, "histograms.jl — equals the sample mean of the raw data; expectations.jl — the g = x quadrature"),
  :var             => (:independent, "histograms.jl — equals the sample variance of the raw data, corrected and uncorrected"),
  :std             => (:independent, "histograms.jl — the square root of the recovered variance against the raw data's"),
  :central_moment  => (:independent, "histograms.jl — equals the direct central moments of the raw data, and the true Gaussian moments"),

  # ---- power-sum quadrature ----
  :expectation                              => (:analytic, "expectations.jl — E[cos X] = e^{-1/2} and the MGF on Gaussian data; polynomial exactness against the raw sum"),
  :expectation_standard_error               => (:independent, "expectations.jl — reduces to s/√n of the raw samples, and the Monte-Carlo spread across independent datasets"),
  :expectation_with_standard_error          => (:analytic, "expectations.jl — value e^{-1/2} and standard error √(((1 + e^{-2})/2 − e^{-1})/n) in closed form on N(0, 1)"),
  :pushforward_density                      => (:analytic, "expectations.jl — the lognormal pushforward of a Gaussian and the non-monotone χ²₁ pushforward, in closed form"),
  :pushforward_density_standard_error       => (:independent, "expectations.jl — the Monte-Carlo spread of the estimate across independent datasets"),
  :pushforward_density_with_standard_error  => (:analytic, "expectations.jl — value the lognormal pdf, standard error the two-term KDE variance at it, in closed form"),

  # ---- the marginal density and its bandwidth ----
  :density                     => (:analytic, "oracles.jl — recovers the Gaussian pdf; density_convergence.jl — the analytic-PDF suite at the theoretical rates"),
  :density_standard_error      => (:independent, "marginaldensity.jl — the Monte-Carlo spread of the estimate across independent datasets, and the oracle coverage of the band"),
  :density_with_standard_error => (:analytic, "marginaldensity.jl — value the N(0, 1 + h²) pdf (the exact KDE mean), standard error the two-term variance at it"),
  :bandwidth                   => (:analytic, "bandwidth.jl, oracles.jl — Silverman's (4/3)^{1/5} σ n^{-1/5} closed form, the Epanechnikov canonical ratio, the degree-2 n^{-1/9} rule"),

  # ---- the smoothing kernels: their constants by an independent quadrature ----
  :gaussian_kernel     => (:analytic, "kernels.jl — ∫K = 1, μ₂(K) = 1, R(K) = 1/(2√π) by an independent quadrature"),
  :epanechnikov_kernel => (:analytic, "kernels.jl — ∫K = 1, μ₂(K) = 1/5, R(K) = 3/5 by an independent quadrature"),

  # ---- conditional expectations (bivariate-Gaussian and heteroscedastic oracles) ----
  :conditional_expectation                              => (:analytic, "conditionalexpectations.jl — E[cos Y | x] of the bivariate Gaussian in closed form"),
  :conditional_expectation_standard_error               => (:independent, "conditionalexpectations.jl — the Monte-Carlo spread of the estimate across independent datasets"),
  :conditional_expectation_with_standard_error          => (:analytic, "conditionalexpectations.jl — value cos(ρx) e^{-v/2}, standard error √(R(K) Var(cos Y|x)/(n h f_X(x))) in closed form"),
  :conditional_mean                                     => (:analytic, "conditionalexpectations.jl — recovers m(x) = ρx of the bivariate Gaussian; conditional_convergence.jl — the n^{-2/5} rate on m(x) = sin x"),
  :conditional_variance                                 => (:analytic, "conditionalexpectations.jl — recovers 1 − ρ² of the bivariate Gaussian and the heteroscedastic σ²(x) oracle"),
  :conditional_variance_standard_error                  => (:independent, "conditionalexpectations.jl — the Monte-Carlo spread across independent datasets, Gaussian and heavy-tailed responses"),
  :conditional_variance_with_standard_error             => (:analytic, "conditionalexpectations.jl — value 1 − ρ², standard error √(ℓ(x)·2σ⁴) with the Gaussian μ₄ = 3σ⁴, in closed form"),
  :conditional_pushforward_density                      => (:analytic, "conditionalexpectations.jl — the lognormal conditional pushforward of the bivariate Gaussian in closed form"),
  :conditional_pushforward_density_standard_error       => (:independent, "conditionalexpectations.jl — the Monte-Carlo spread across independent datasets and the oracle coverage of the band"),
  :conditional_pushforward_density_with_standard_error  => (:analytic, "conditionalexpectations.jl — value the conditional lognormal pdf, standard error the two-term variance at it, in closed form"),

  # ---- conditional and transition densities (bivariate-Gaussian and Ornstein–Uhlenbeck oracles) ----
  :conditional_density                     => (:analytic, "conditionaldensity.jl — the bivariate-Gaussian conditional N(ρx, 1 − ρ²); conditional_convergence.jl — the n^{-2/3} rate"),
  :conditional_density_standard_error      => (:independent, "conditionaldensity.jl — the Monte-Carlo spread of the estimate across independent datasets"),
  :conditional_density_with_standard_error => (:analytic, "conditionaldensity.jl — value N(ρx, 1 − ρ²), standard error the two-term double-kernel variance at it, in closed form"),
  :transition_density                      => (:analytic, "oracles.jl — the Ornstein–Uhlenbeck transition density N(x e^{-θτ}, s²(1 − e^{-2θτ})) in closed form"),
  :transition_density_standard_error       => (:analytic, "oracles.jl — the two-term double-kernel variance with the OU transition law and stationary marginal, in closed form"),
  :transition_density_with_standard_error  => (:analytic, "oracles.jl — value the OU transition law, standard error its two-term variance, in closed form"),

  # ---- HDF5 persistence (the weak-dependency extension) ----
  :write_histogram => (:independent, "hdf5.jl — the group's datasets and attributes are read back with raw HDF5 calls and compared to the summary's arrays and parameters"),
  :read_histogram  => (:independent, "hdf5.jl — the reload is bit-exact against the written summary's arrays, and the result still merges"),
)

GBGCodingTools.TestHarness.oracle_registry(GBGStreamingStatistics, PACKAGE_MODULES,
    Dict("operation" => [Dict("name" => String(op), "class" => String(class), "evidence" => evidence)
                         for (op, (class, evidence)) in ORACLE_PROVENANCE]),
    @__DIR__)
