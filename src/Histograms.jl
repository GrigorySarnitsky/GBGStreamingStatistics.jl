# Histograms — the carriers (universal lower half): the finely binned empirical measures that are the
# canonical summary `T`. Linear (WARPing) binning, per-bin counts, and per-bin central moments via the
# Pébay recursion. Rationale: `docs/ARCHITECTURE.md` §3.

module Histograms

import ..Accumulators: observe!
import Base: merge!
using ..Accumulators: AbstractAccumulator, _observe_moment!, _merge_moment!

export Histogram1D, Histogram2D, bin_index, observe_pairs!


"""
    Histogram1D <: AbstractAccumulator
    Histogram1D(lo, hi, nbins; order = 4, accumulator = Float64)

A **finely binned empirical measure** of a conditioning variable ``X`` — the canonical summary
``T(\\hat\\mu_n)`` (`docs/ARCHITECTURE.md` §3). `nbins` equal-width (WARPing) bins tile ``[lo, hi)``;
choosing the range to cover the data (e.g. ``\\pm 5\\sigma``) with ``\\gtrsim 10^3`` bins makes the
bin width ``b \\approx \\sigma/100 \\ll h`` so the summary is *approximately sufficient* (the
``O(b^2)`` binning error is dominated by the ``O(h^2)`` smoothing bias).

Each bin carries its **count**, a running **mean**, and the unnormalized **central moments**
``M_p = \\sum_{i \\in \\text{bin}}(X_i - \\bar X)^p`` for ``p = 2 \\dots`` `order`, accumulated by the
Pébay recursion. Counts are `UInt64`; the mean and moments are accumulated in `accumulator`
(`Float64` by default) regardless of the input element type — *accumulate wider than the data*. The
default `order = 4` reaches kurtosis; the heavy-tailed channels (e.g. the Pawula ``Q``-criterion's
``M_6``, whose sampling variance reaches ``M_{12}``) construct with a higher `order`.

Samples outside ``[lo, hi)`` are currently **skipped** (not yet counted) — set the range to cover the
data; out-of-range bookkeeping is a recorded gap (`CLAUDE.md` roadmap).
"""
struct Histogram1D{T<:AbstractFloat,C<:AbstractVector,V<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  lo::T          # lower edge of the binned range
  hi::T          # upper edge of the binned range
  scale::T       # nbins / (hi − lo), precomputed: bin of x is floor((x − lo)·scale) + 1
  counts::C      # length nbins, per-bin integer counts (UInt64)
  means::V       # length nbins, per-bin running mean (accumulator eltype)
  moments::M     # (order−1) × nbins, per-bin central moments M_2 … M_order (accumulator eltype)
end


function Histogram1D(lo::Real, hi::Real, nbins::Integer;
                     order::Integer = 4, accumulator::Type{<:AbstractFloat} = Float64)
  hi > lo || throw(ArgumentError("Histogram1D requires hi > lo (got lo = $lo, hi = $hi)"))
  nbins ≥ 1 || throw(ArgumentError("Histogram1D requires nbins ≥ 1 (got $nbins)"))
  order ≥ 2 || throw(ArgumentError("Histogram1D requires order ≥ 2 (got $order)"))
  loT, hiT = promote(float(lo), float(hi))
  scale = oftype(loT, nbins) / (hiT - loT)
  counts = zeros(UInt64, nbins)
  means = zeros(accumulator, nbins)
  moments = zeros(accumulator, order - 1, nbins)
  return Histogram1D(loT, hiT, scale, counts, means, moments)
end


_nbins(h::Histogram1D) = length(h.counts)
_order(h::Histogram1D) = size(h.moments, 1) + 1
_bincenter(h::Histogram1D, i::Integer) = h.lo + (oftype(h.lo, i) - oftype(h.lo, 0.5)) / h.scale
_total(h::Histogram1D) = sum(h.counts)


"""
    bin_index(summary, x) -> i

The **binning map**: the index ``i`` of the bin containing the value `x`. For a
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) over ``[lo, hi)`` with `nbins`
bins, ``i = \\lfloor (x - lo)\\,\\mathrm{scale} \\rfloor + 1``; a valid bin satisfies
``1 \\le i \\le`` `nbins`, and `x` outside the range yields ``i < 1`` or ``i >`` `nbins`. Exposed so a
consumer can pre-quantize a snapshot buffer (e.g. store compact `UInt16` bin indices) before forming
lagged pairs for [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!).
"""
@inline bin_index(h::Histogram1D, x) = floor(Int, (x - h.lo) * h.scale) + 1


function observe!(h::Histogram1D, samples)
  nb = _nbins(h)
  ord = _order(h)
  one_count = one(eltype(h.counts))
  @inbounds for x in samples
    i = floor(Int, (x - h.lo) * h.scale) + 1
    (1 ≤ i ≤ nb) || continue
    n0 = h.counts[i]
    _observe_moment!(h.means, h.moments, i, n0, x, ord)
    h.counts[i] = n0 + one_count
  end
  return h
end


function observe!(h::Histogram1D, x::Real)
  i = floor(Int, (x - h.lo) * h.scale) + 1
  (1 ≤ i ≤ _nbins(h)) || return h
  @inbounds n0 = h.counts[i]
  _observe_moment!(h.means, h.moments, i, n0, x, _order(h))
  @inbounds h.counts[i] = n0 + one(eltype(h.counts))
  return h
end


_compatible(a::Histogram1D, b::Histogram1D) =
  a.lo == b.lo && a.hi == b.hi && _nbins(a) == _nbins(b) && _order(a) == _order(b)


function merge!(into::Histogram1D, from::Histogram1D)
  _compatible(into, from) ||
    throw(ArgumentError("merge!: histograms must share lo, hi, bin count, and moment order"))
  ord = _order(into)
  @inbounds for b in 1:_nbins(into)
    nA = into.counts[b]
    nB = from.counts[b]
    _merge_moment!(into.means, into.moments, b, nA, from.means, from.moments, b, nB, ord)
    into.counts[b] = nA + nB
  end
  return into
end


# Global (over all bins) central moments, by folding the per-bin moment-sets through the same Pébay
# merge (with the running accumulator held in a one-element carrier). Returns (n, mean, M) where
# `M[p-1] = Σ_i (x_i − x̄)^p`. Internal — used by bandwidth selection and the test oracles.
function _global_moments(h::Histogram1D)
  ord = _order(h)
  A = eltype(h.means)
  gmean = zeros(A, 1)
  gmom = zeros(A, ord - 1, 1)
  ng = 0.0
  @inbounds for b in 1:_nbins(h)
    nB = h.counts[b]
    _merge_moment!(gmean, gmom, 1, ng, h.means, h.moments, b, nB, ord)
    ng += nB
  end
  return ng, gmean[1], gmom
end

_sample_mean(h::Histogram1D) = _global_moments(h)[2]

function _sample_central_moment(h::Histogram1D, p::Integer)
  2 ≤ p ≤ _order(h) ||
    throw(ArgumentError("central-moment order p must satisfy 2 ≤ p ≤ $(_order(h)) (got $p)"))
  ng, _, gmom = _global_moments(h)
  return iszero(ng) ? oftype(gmom[1], NaN) : gmom[p - 1, 1] / ng
end

function _sample_variance(h::Histogram1D)
  ng, _, gmom = _global_moments(h)
  return iszero(ng) ? oftype(gmom[1], NaN) : gmom[1, 1] / ng
end

_sample_std(h::Histogram1D) = sqrt(_sample_variance(h))


"""
    Histogram2D <: AbstractAccumulator

A **joint binned empirical measure** of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities (`docs/ARCHITECTURE.md` §3, §4.3), with the same per-bin carriers as
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) on the product grid. For a
two-time quantity the conditioning value ``X(t)`` and the response ``Y(t+\\tau)`` must be paired *per
source*: the entry is made at the response step into the bin of ``X(t)`` (two separate marginals
cannot recover the joint — for the Ornstein–Uhlenbeck oracle both marginals are the same stationary
Gaussian and the entire ``\\tau``-dependence lives in the pair). The lag ``\\tau`` is an explicit
label on the carrier.

!!! note "Forthcoming"
    Declared; the two-time carrier and [`observe_pairs!`](@ref
    GBGStreamingStatistics.Histograms.observe_pairs!) land in the next pass, after the `Histogram1D`
    foundation (`CLAUDE.md` roadmap).
"""
struct Histogram2D{T,C<:AbstractArray,A<:AbstractArray} <: AbstractAccumulator
  xlo::T         # lower edge of the X range
  xhi::T         # upper edge of the X range
  ylo::T         # lower edge of the Y range
  yhi::T         # upper edge of the Y range
  lag::Int       # the lag τ this joint carrier is accumulated at (0 for the static joint)
  counts::C      # per-cell integer counts on the (X, Y) product grid
  moments::A     # per-cell central moments
end


"""
    observe_pairs!(summary, x_past, y_now; lag) -> summary

Fold a batch of **pre-formed, lag-labelled pairs** ``\\bigl(X(t),\\, Y(t+\\tau)\\bigr)`` into the
per-lag carrier in place — the memoryless surface for two-time statistics (conditional and transition
densities, conditional increment moments; `docs/ARCHITECTURE.md` §3, §4.3–4.4). Everything that
*produces* the pairs — the snapshot buffer of past ``X(t)``, which sources to sample, which lags to
form, the origin cadence — is the **consumer's** responsibility; this package only bins the pairs it
is handed, keeping the small per-lag histograms and never the ``O(ML)`` snapshot buffer.

!!! note "Forthcoming"
    Declared; lands with [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) in the
    next pass.
"""
function observe_pairs! end

end # module Histograms
