# Accumulators — the foundation (universal lower half): the mergeable-monoid accumulator interface
# and the per-bin central-moment recursions. Pure numerics over a stream; transfers verbatim to every
# consumer. The conceptual model and conventions live in `docs/ARCHITECTURE.md` (rationale) and
# `CLAUDE.md` (operational rules). Rationale: §3, §5.

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
additive monoid over the stream; see [`merge!`](@ref Base.merge!)), and **rich enough** that every
reported estimator is recoverable from it post hoc, with no second pass over the data. The summary is
accumulated online and is *response-agnostic*; all smoothing (kernel, bandwidth, channel) is applied
post hoc. Concrete carriers ([`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D),
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

Allocation-free on a preallocated `summary` (required for GPU hot loops, where a per-call device
temporary escapes the host GC and wedges the run). Per-bin central moments are accumulated by the
**Pébay (2008) recursion** in `Float64` regardless of input element type — *accumulate wider than the
data* (`docs/ARCHITECTURE.md` §3). The scatter-add binning is a serial CPU loop here; the
GPU-resident KernelAbstractions + Atomix kernel (per-block privatization, merged via the monoid) and
Kahan–Neumaier compensation for the very-large-``n`` regime are forthcoming (`CLAUDE.md` roadmap).

Reference: P. Pébay (2008), *Formulas for robust, one-pass parallel computation of covariances and
arbitrary-order statistical moments*, Sandia Report SAND2008-6212.
"""
function observe! end


"""
    merge!(into::AbstractAccumulator, from::AbstractAccumulator) -> into

Combine two summaries into `into`, forming the summary of the concatenated stream — the monoid
operation that makes accumulation additive: ``T(\\hat\\mu_{A \\cup B}) = T(\\hat\\mu_A) \\oplus
T(\\hat\\mu_B)``. Counts add; the per-bin central moments combine by the **Pébay parallel-merge
formula** (of which the single-sample [`observe!`](@ref
GBGStreamingStatistics.Accumulators.observe!) update is the ``n_B = 1`` case, so the two agree
exactly). This is what makes multi-device reduction and time-chunk stitching free, and lets the hot
loop privatize (per-block local histograms, merged hierarchically).

`merge!` is intrinsically a finite-sample operation, so it has no population estimand; its
specification is the monoid law, tested as *accumulate-then-merge `=` accumulate-together* (to
round-off). Extends `Base.merge!`.
"""
merge!(into::AbstractAccumulator, from::AbstractAccumulator) =
  error("merge! is not implemented for $(typeof(into))")


# ---------------------------------------------------------------------------------------------------
# Pébay per-bin central-moment recursions — the hard, write-once numerics (`docs/ARCHITECTURE.md` §3).
#
# A carrier stores, per bin `b`: a count, a running mean `means[b]`, and the unnormalized central
# moments `moments[p-1, b] = M_p = Σ_{i ∈ bin}(x_i − mean)^p` for `p = 2 … order`. Both updaters write
# the moments high-order → low-order so each `M_p` reads the *pre-update* lower moments it needs;
# counts are owned by the carrier and updated by the caller.
# ---------------------------------------------------------------------------------------------------

"""
    _observe_moment!(means, moments, b, n0, x, order)

Fold the value `x` into bin `b`, given the bin's prior count `n0`. Updates `means[b]` and
`moments[2..order, b]` in place (the caller then increments the count). This is the ``n_B = 1`` case
of [`_merge_moment!`](@ref), so streaming accumulation and `merge!` produce identical moments.
"""
@inline function _observe_moment!(means, moments, b::Integer, n0::Real, x, order::Integer)
  if iszero(n0)
    @inbounds means[b] = x
    return nothing
  end
  n0f = float(n0)
  n1 = n0f + one(n0f)
  @inbounds μ = means[b]
  δ = x - μ
  δn = δ / n1
  for p in order:-1:2
    @inbounds Mp = moments[p - 1, b]
    acc = zero(Mp)
    term = one(δn)                                   # accumulates (−δn)^k
    for k in 1:(p - 2)
      term *= -δn
      @inbounds acc += binomial(p, k) * moments[p - k - 1, b] * term
    end
    acc += δ^p * (n0f / n1)^p * (one(n0f) - (-one(n0f) / n0f)^(p - 1))
    @inbounds moments[p - 1, b] = Mp + acc
  end
  @inbounds means[b] = μ + δn
  return nothing
end


"""
    _merge_moment!(means_d, moments_d, bd, nA, means_s, moments_s, bs, nB, order)

Merge source bin `bs` (count `nB`) into destination bin `bd` (count `nA`) by the Pébay parallel-merge
formula, updating `means_d[bd]` and `moments_d[2..order, bd]` in place to represent the union (the
caller updates the count). `moments_s` is read-only, so the same source may be merged into many
destinations (used by the global-moment reduction with `bd = 1`).
"""
@inline function _merge_moment!(means_d, moments_d, bd::Integer, nA::Real,
                                means_s, moments_s, bs::Integer, nB::Real, order::Integer)
  iszero(nB) && return nothing
  if iszero(nA)
    @inbounds means_d[bd] = means_s[bs]
    for j in 1:(order - 1)
      @inbounds moments_d[j, bd] = moments_s[j, bs]
    end
    return nothing
  end
  nAf = float(nA); nBf = float(nB); n = nAf + nBf
  @inbounds δ = means_s[bs] - means_d[bd]
  a = -nBf * δ / n
  c = nAf * δ / n
  for p in order:-1:2
    @inbounds acc = moments_d[p - 1, bd] + moments_s[p - 1, bs]
    ak = one(a); ck = one(c)
    for k in 1:(p - 2)
      ak *= a; ck *= c
      bk = binomial(p, k)
      @inbounds acc += bk * (moments_d[p - k - 1, bd] * ak + moments_s[p - k - 1, bs] * ck)
    end
    acc += δ^p * (nAf * nBf / n)^p * (one(nBf) / nBf^(p - 1) - (-one(nAf) / nAf)^(p - 1))
    @inbounds moments_d[p - 1, bd] = acc
  end
  @inbounds means_d[bd] += δ * nBf / n
  return nothing
end

end # module Accumulators
