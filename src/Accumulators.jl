# Accumulators — the foundation (universal lower half): the mergeable-monoid accumulator interface
# and the backend-agnostic accumulation primitives. Pure numerics over a stream; transfers verbatim
# to every consumer. The conceptual model and conventions are not duplicated here — they live in
# `docs/ARCHITECTURE.md` (the rationale) and `CLAUDE.md` (the operational rules). Rationale: §3, §5.

module Accumulators

using KernelAbstractions
using Atomix

import Base: merge!

export AbstractAccumulator, observe!, merge!


"""
    AbstractAccumulator

Supertype of every **mergeable summary** ``T(\\hat\\mu_n)`` of the empirical measure
``\\hat\\mu_n = \\frac1n \\sum_{i=1}^{n} \\delta_{X_i}`` of an observed stream ``X_1, X_2, \\dots``
(``\\delta_x`` the unit point mass at ``x``, ``n`` the number of samples seen).

A summary is **finite-dimensional** — the ``n`` point masses are never stored — **mergeable** (an
additive monoid over the stream; see [`merge!`](@ref Base.merge!)), and **rich enough** that every reported
estimator is recoverable from it post hoc, with no second pass over the data. The summary is
accumulated online and is *response-agnostic*; all smoothing (kernel, bandwidth, channel) is applied
post hoc. Concrete carriers ([`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D), [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D)) subtype this.

!!! note "Scaffold stage"
    Declared; the accumulation numerics are not yet implemented.
"""
abstract type AbstractAccumulator end


"""
    observe!(summary::AbstractAccumulator, batch) -> summary

Fold a `batch` of samples into `summary` in place, advancing the accumulated summary
``T(\\hat\\mu_n)`` of the empirical measure (see [`AbstractAccumulator`](@ref GBGStreamingStatistics.Accumulators.AbstractAccumulator)). Accumulation is
**response-agnostic**: it knows nothing of bandwidth, kernel, or which estimator will later be asked
for — the architectural invariant that delivers *no reruns*.

Allocation-free on a preallocated `summary` (required for GPU hot loops, where a per-call device
temporary escapes the host GC and wedges the run) and backend-agnostic: the atomic scatter-add
binning is one KernelAbstractions + Atomix kernel that runs on CPU and any GPU backend. Per-bin
central moments are accumulated by the Pébay recursion with Kahan–Neumaier compensated summation, in
`Float64` regardless of input element type — *accumulate wider than the data* (`docs/ARCHITECTURE.md`
§3).

Reference: P. Pébay (2008), *Formulas for robust, one-pass parallel computation of covariances and
arbitrary-order statistical moments*, Sandia Report SAND2008-6212.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function observe! end


"""
    merge!(into::AbstractAccumulator, from::AbstractAccumulator) -> into

Combine two summaries by adding their **unnormalized** contents — counts and power-sums add — the
monoid operation that makes accumulation additive over the stream:
``T(\\hat\\mu_{A \\cup B}) = T(\\hat\\mu_A) \\oplus T(\\hat\\mu_B)``. This is just additivity of
measures, and it is what makes multi-device reduction and time-chunk stitching free, and what lets
the hot loop privatize (per-block local histograms, merged hierarchically).

`merge!` is intrinsically a finite-sample operation, so it has no population estimand; its
specification is the monoid law, tested as *accumulate-then-merge `=` accumulate-together* (to
round-off). Extends `Base.merge!`.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
merge!(into::AbstractAccumulator, from::AbstractAccumulator) =
  error("merge! is not yet implemented (scaffold stage)")


# The backend-agnostic scatter-add primitive (one KernelAbstractions + Atomix kernel) and the
# Pébay + Kahan per-bin central-moment recursion live here once implemented — the hard, write-once
# numerics the whole package depends on (`docs/ARCHITECTURE.md` §3, §5).

end # module Accumulators
