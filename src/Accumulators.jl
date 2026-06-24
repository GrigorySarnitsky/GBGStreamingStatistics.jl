# Accumulators — the foundation (universal lower half): the mergeable-monoid accumulator interface.
# The concrete summaries and their KernelAbstractions + Atomix scatter kernels live in `Histograms`.
# Conceptual model and conventions: `docs/ARCHITECTURE.md` (rationale), `CLAUDE.md` (operational rules).

module Accumulators

import Base: merge!

export AbstractAccumulator, observe!, merge!


"""
    AbstractAccumulator

Supertype of every **mergeable summary** ``T(\\hat\\mu_n)`` of the [empirical measure](https://en.wikipedia.org/wiki/Empirical_measure)
``\\hat\\mu_n = \\frac1n \\sum_{i=1}^{n} \\delta_{X_i}`` of an observed stream ``X_1, X_2, \\dots``
(``\\delta_x`` the unit point mass at ``x``, ``n`` the number of samples seen).

A summary is **finite-dimensional** — the ``n`` point masses are never stored — **mergeable** (an
additive monoid over the stream; see [`merge!`](@ref Base.merge!)), and **rich enough** that every
reported estimator is recoverable from it post hoc, with no second pass over the data. The summary is
accumulated online and is *response-agnostic*; all smoothing (kernel, bandwidth, channel) is applied
post hoc. Concrete summaries ([`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D),
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D)) subtype this.
"""
abstract type AbstractAccumulator end


"""
    observe!(summary::AbstractAccumulator, batch) -> summary

Fold a `batch` of samples into `summary` in place, advancing the accumulated summary
``T(\\hat\\mu_n)`` of the empirical measure (see
[`AbstractAccumulator`](@ref GBGStreamingStatistics.Accumulators.AbstractAccumulator)). Accumulation
is **response-agnostic**: it knows nothing of bandwidth, kernel, or which estimator will later be
asked for — the architectural invariant that delivers *no reruns*.

The binning is an **atomic scatter-add**: one KernelAbstractions + Atomix kernel — one work-item per
sample — accumulates, per bin, the **count** and the **power-sums**
``S_p = \\sum_{i\\in\\text{bin}}(X_i - c_j)^p`` about the bin center ``c_j``, directly into the
preallocated summary (no scratch). It runs identically on a CPU `Array` and any GPU backend, and is
allocation-free up to the kernel launch.
Because ``X_i - c_j`` is small within a bin, the local power-sums are well-conditioned (no compensated
summation needed); global central moments are recovered post hoc by the shift formula
(`docs/ARCHITECTURE.md` §3). Accumulated in the summary's accumulator type (`Float64` by default),
never narrowed to the input element type — *accumulate wider than the data*.

Reference: D. W. Scott (2015), *Multivariate Density Estimation*, 2nd ed.; M. D. Cattaneo et al.
(2024), *On binscatter* (the local-moment / WARPing construction).
"""
function observe! end


"""
    merge!(into::AbstractAccumulator, from::AbstractAccumulator) -> into

Combine two summaries into `into`, forming the summary of the concatenated stream — the monoid
operation that makes accumulation additive: ``T(\\hat\\mu_{A \\cup B}) = T(\\hat\\mu_A) \\oplus
T(\\hat\\mu_B)``. Because compatible summaries share the same binning (hence the same bin centers),
their **counts and power-sums simply add elementwise**, so `merge!` is exact and trivially parallel —
this is what makes multi-device reduction and time-chunk stitching free, and lets the hot loop
privatize (per-block local histograms, merged hierarchically).

`merge!` is intrinsically a finite-sample operation, so it has no population estimand; its
specification is the monoid law, tested as *accumulate-then-merge `=` accumulate-together*. Extends
`Base.merge!`.
"""
merge!(into::AbstractAccumulator, from::AbstractAccumulator) =
  error("merge! is not implemented for $(typeof(into))")

end # module Accumulators
