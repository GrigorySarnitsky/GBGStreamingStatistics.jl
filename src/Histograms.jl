# Histograms — the summaries (universal lower half): the finely binned empirical measures that are the
# canonical summary `T`. Linear (WARPing) binning; per-bin counts + power-sums about the bin center,
# accumulated by one portable KernelAbstractions + Atomix scatter kernel. Rationale: `docs/ARCHITECTURE.md` §3.

module Histograms

import ..Accumulators: observe!
import Base: merge!
using ..Accumulators: AbstractAccumulator
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
import Atomix

export Histogram1D, Histogram2D, bin_index, observe_pairs!


"""
    Histogram1D <: AbstractAccumulator
    Histogram1D(lo, hi, nbins; order = 4, accumulator = Float64)

A **finely binned empirical measure** of a conditioning variable ``X`` — the canonical summary
``T(\\hat\\mu_n)`` (`docs/ARCHITECTURE.md` §3). `nbins` equal-width (WARPing) bins tile ``[lo, hi)``;
choosing the range to cover the data (e.g. ``\\pm 5\\sigma``) with ``\\gtrsim 10^3`` bins makes the
bin width ``b \\approx \\sigma/100 \\ll h`` so the summary is *approximately sufficient*.

Each bin carries its **count** and the **power-sums** ``S_p = \\sum_{i \\in \\text{bin}}(X_i - c_j)^p``
about the fixed bin center ``c_j``, for ``p = 1 \\dots`` `order`. Because ``X_i - c_j`` is small within
a bin these are well-conditioned, so accumulation is a plain atomic scatter-add (`observe!`) that runs
identically on CPU and GPU, and `merge!` is elementwise addition. The power-sums use the *actual*
``X_i`` (not the center), so the global central moments recovered post hoc (by the shift formula) are
exact; they are stable for the even high-order moments (the Pawula ``Q``-criterion's ``M_4``, ``M_6``).
Counts are `UInt64`; power-sums are accumulated in `accumulator` (`Float64` by default) regardless of
the input element type — *accumulate wider than the data*.

Samples outside ``[lo, hi)`` are **skipped** (not counted) — set the range to cover the data.
"""
struct Histogram1D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  lo::T          # lower edge of the binned range
  hi::T          # upper edge of the binned range
  scale::T       # nbins / (hi − lo); bin of x is floor((x − lo)·scale) + 1
  counts::C      # length nbins, per-bin integer counts (UInt64)
  powersums::M   # order × nbins, per-bin power-sums S_1 … S_order about the bin center
end


function Histogram1D(lo::Real, hi::Real, nbins::Integer;
                     order::Integer = 4, accumulator::Type{<:AbstractFloat} = Float64)
  hi > lo || throw(ArgumentError("Histogram1D requires hi > lo (got lo = $lo, hi = $hi)"))
  nbins ≥ 1 || throw(ArgumentError("Histogram1D requires nbins ≥ 1 (got $nbins)"))
  order ≥ 2 || throw(ArgumentError("Histogram1D requires order ≥ 2 (got $order)"))
  loT, hiT = promote(float(lo), float(hi))
  scale = oftype(loT, nbins) / (hiT - loT)
  counts = zeros(UInt64, nbins)
  powersums = zeros(accumulator, order, nbins)
  return Histogram1D(loT, hiT, scale, counts, powersums)
end


_nbins(h::Histogram1D) = length(h.counts)
_order(h::Histogram1D) = size(h.powersums, 1)
_bincenter(h::Histogram1D, i::Integer) = h.lo + (oftype(h.lo, i) - oftype(h.lo, 0.5)) / h.scale
_total(h::Histogram1D) = sum(h.counts)


"""
    bin_index(summary, x) -> i

The **binning map**: the index ``i`` of the bin containing `x`. For a
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) over ``[lo, hi)`` with `nbins`
bins, ``i = \\lfloor (x - lo)\\,\\mathrm{scale} \\rfloor + 1``; a valid bin satisfies
``1 \\le i \\le`` `nbins`, and `x` outside the range yields ``i < 1`` or ``i >`` `nbins`. Exposed so a
consumer can pre-quantize a snapshot buffer (e.g. store compact `UInt16` bin indices) before forming
lagged pairs for [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!).
"""
@inline bin_index(h::Histogram1D, x) = floor(Int, (x - h.lo) * h.scale) + 1


# The portable scatter-add: one work-item per sample atomically accumulates the count and the
# per-bin power-sums about the bin center, directly into the preallocated summary. Runs on CPU and
# any GPU backend (KernelAbstractions + Atomix). Out-of-range samples are skipped.
@kernel function _scatter1d!(counts, powersums, @Const(samples), lo, scale, nbins, order)
  s = @index(Global, Linear)
  if s ≤ length(samples)
    @inbounds begin
      x = samples[s]
      i = floor(Int, (x - lo) * scale) + 1
      if 1 ≤ i ≤ nbins
        Atomix.@atomic counts[i] += one(eltype(counts))
        center = lo + (oftype(lo, i) - oftype(lo, 0.5)) / scale
        dx = x - center
        base = (i - 1) * order
        dxp = dx
        for p in 1:order
          Atomix.@atomic powersums[base + p] += dxp
          dxp = dxp * dx
        end
      end
    end
  end
end


function observe!(h::Histogram1D, samples)
  backend = get_backend(h.counts)
  _scatter1d!(backend)(h.counts, h.powersums, samples, h.lo, h.scale, _nbins(h), _order(h);
                       ndrange = length(samples))
  synchronize(backend)
  return h
end


# Scalar fold — a direct host update (scalar-indexes the summary, so host arrays only; batch on GPU).
function observe!(h::Histogram1D, x::Real)
  i = floor(Int, (x - h.lo) * h.scale) + 1
  (1 ≤ i ≤ _nbins(h)) || return h
  @inbounds h.counts[i] += one(eltype(h.counts))
  dx = x - _bincenter(h, i)
  dxp = dx
  @inbounds for p in 1:_order(h)
    h.powersums[p, i] += dxp
    dxp *= dx
  end
  return h
end


_compatible(a::Histogram1D, b::Histogram1D) =
  a.lo == b.lo && a.hi == b.hi && _nbins(a) == _nbins(b) && _order(a) == _order(b)


function merge!(into::Histogram1D, from::Histogram1D)
  _compatible(into, from) ||
    throw(ArgumentError("merge!: histograms must share lo, hi, bin count, and order"))
  into.counts .+= from.counts          # counts and power-sums add (same binning ⇒ same centers)
  into.powersums .+= from.powersums
  return into
end


# Global (over all bins) central moments, recovered host-side from the per-bin power-sums by the
# shift formula. Returns (n, mean, M) with `M[p-1] = Σ_i (x_i − x̄)^p` for p = 2 … order. Reads the
# summaries to the host first (a no-op for CPU arrays). Internal — used by bandwidth and the oracles.
function _global_moments(h::Histogram1D)
  P = _order(h)
  counts = Array(h.counts)
  S = Array(h.powersums)
  nbins = length(counts)
  n = 0.0
  sx = 0.0
  for j in 1:nbins
    nj = counts[j]
    iszero(nj) && continue
    n += nj
    sx += _bincenter(h, j) * nj + S[1, j]      # Σ x = Σ (c_j n_j + S_1[j])
  end
  M = zeros(P - 1)
  iszero(n) && return 0.0, NaN, fill!(M, NaN)
  xbar = sx / n
  for j in 1:nbins
    nj = counts[j]
    iszero(nj) && continue
    d = _bincenter(h, j) - xbar
    for p in 2:P
      acc = 0.0
      for k in 0:p
        Sk = k == 0 ? float(nj) : S[k, j]
        acc += binomial(p, k) * d^(p - k) * Sk
      end
      M[p - 1] += acc
    end
  end
  return n, xbar, M
end

_sample_mean(h::Histogram1D) = _global_moments(h)[2]

function _sample_central_moment(h::Histogram1D, p::Integer)
  2 ≤ p ≤ _order(h) ||
    throw(ArgumentError("central-moment order p must satisfy 2 ≤ p ≤ $(_order(h)) (got $p)"))
  n, _, M = _global_moments(h)
  return iszero(n) ? oftype(M[1], NaN) : M[p - 1] / n
end

_sample_variance(h::Histogram1D) = _sample_central_moment(h, 2)
_sample_std(h::Histogram1D) = sqrt(_sample_variance(h))


"""
    Histogram2D <: AbstractAccumulator
    Histogram2D(xlo, xhi, nx, ylo, yhi, ny; lag = 0, accumulator = Float64)

A **joint binned empirical measure** of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities (`docs/ARCHITECTURE.md` §3, §4.3). `nx × ny` equal-width cells tile
``[xlo, xhi) \\times [ylo, yhi)``; each cell `(i_x, i_y)` (flattened `i_x + (i_y-1)·nx`) carries its
**count** and the order-2 bivariate **power-sums** ``S_{rs} = \\sum_{(X,Y)\\in\\text{cell}} (X -
c^x)^r (Y - c^y)^s`` about the cell center ``(c^x, c^y)`` — ``S_{10}, S_{01}, S_{20}, S_{11},
S_{02}`` — the local sufficient statistic for the joint density and for local-linear conditioning.
For a two-time quantity ``X(t)`` and the response ``Y(t+\\tau)`` are paired *per source*; the lag
``\\tau`` is an explicit label on the summary. Pairs outside the box are skipped. Same atomic-scatter
accumulation, `merge!`-by-addition, and `Float64` accumulation as
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D).
"""
struct Histogram2D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  xlo::T; xhi::T; xscale::T     # X range and nx/(xhi−xlo)
  ylo::T; yhi::T; yscale::T     # Y range and ny/(yhi−ylo)
  nx::Int
  ny::Int
  lag::Int                       # the lag τ this joint summary is accumulated at (0 for the static joint)
  counts::C                      # length nx·ny, per-cell integer counts (UInt64)
  powersums::M                   # 5 × (nx·ny): rows S_10, S_01, S_20, S_11, S_02 about the cell center
end


function Histogram2D(xlo::Real, xhi::Real, nx::Integer, ylo::Real, yhi::Real, ny::Integer;
                     lag::Integer = 0, accumulator::Type{<:AbstractFloat} = Float64)
  (xhi > xlo && yhi > ylo) ||
    throw(ArgumentError("Histogram2D requires xhi > xlo and yhi > ylo"))
  (nx ≥ 1 && ny ≥ 1) || throw(ArgumentError("Histogram2D requires nx ≥ 1 and ny ≥ 1"))
  xloT, xhiT, yloT, yhiT = promote(float(xlo), float(xhi), float(ylo), float(yhi))
  T = typeof(xloT)
  xscale = T(nx) / (xhiT - xloT)
  yscale = T(ny) / (yhiT - yloT)
  counts = zeros(UInt64, nx * ny)
  powersums = zeros(accumulator, 5, nx * ny)
  return Histogram2D(xloT, xhiT, xscale, yloT, yhiT, yscale, Int(nx), Int(ny), Int(lag), counts, powersums)
end


_ncells(h::Histogram2D) = length(h.counts)
_total(h::Histogram2D) = sum(h.counts)
_xbincenter(h::Histogram2D, ix::Integer) = h.xlo + (oftype(h.xlo, ix) - oftype(h.xlo, 0.5)) / h.xscale
_ybincenter(h::Histogram2D, iy::Integer) = h.ylo + (oftype(h.ylo, iy) - oftype(h.ylo, 0.5)) / h.yscale


# The 2D scatter-add: one work-item per pair atomically accumulates the per-cell count and the five
# order-2 bivariate power-sums about the cell center. Out-of-box pairs are skipped.
@kernel function _scatter2d!(counts, powersums, @Const(xs), @Const(ys),
                             xlo, xscale, nx, ylo, yscale, ny)
  s = @index(Global, Linear)
  if s ≤ length(xs)
    @inbounds begin
      x = xs[s]
      y = ys[s]
      ix = floor(Int, (x - xlo) * xscale) + 1
      iy = floor(Int, (y - ylo) * yscale) + 1
      if 1 ≤ ix ≤ nx && 1 ≤ iy ≤ ny
        cell = ix + (iy - 1) * nx
        Atomix.@atomic counts[cell] += one(eltype(counts))
        dx = x - (xlo + (oftype(xlo, ix) - oftype(xlo, 0.5)) / xscale)
        dy = y - (ylo + (oftype(ylo, iy) - oftype(ylo, 0.5)) / yscale)
        base = (cell - 1) * 5
        Atomix.@atomic powersums[base + 1] += dx
        Atomix.@atomic powersums[base + 2] += dy
        Atomix.@atomic powersums[base + 3] += dx * dx
        Atomix.@atomic powersums[base + 4] += dx * dy
        Atomix.@atomic powersums[base + 5] += dy * dy
      end
    end
  end
end


"""
    observe_pairs!(summary::Histogram2D, x_past, y_now) -> summary

Fold a batch of **pre-formed pairs** ``\\bigl(X(t),\\, Y(t+\\tau)\\bigr)`` into the per-lag summary in
place (`x_past` ``= X(t)``, `y_now` ``= Y(t+\\tau)``, equal length) — the memoryless surface for
two-time statistics (`docs/ARCHITECTURE.md` §3, §4.3–4.4). The pairs are binned by an atomic
scatter-add (one work-item per pair, CPU and GPU alike); everything that *produces* the pairs — the
snapshot buffer of past ``X(t)``, which sources to sample, which lags to form, the origin cadence — is
the **consumer's** responsibility, formed at the summary's lag ``\\tau``.
"""
function observe_pairs!(h::Histogram2D, x_past, y_now)
  length(x_past) == length(y_now) ||
    throw(ArgumentError("observe_pairs!: x_past and y_now must have equal length"))
  backend = get_backend(h.counts)
  _scatter2d!(backend)(h.counts, h.powersums, x_past, y_now,
                       h.xlo, h.xscale, h.nx, h.ylo, h.yscale, h.ny; ndrange = length(x_past))
  synchronize(backend)
  return h
end


_compatible(a::Histogram2D, b::Histogram2D) =
  a.xlo == b.xlo && a.xhi == b.xhi && a.ylo == b.ylo && a.yhi == b.yhi &&
  a.nx == b.nx && a.ny == b.ny && a.lag == b.lag

function merge!(into::Histogram2D, from::Histogram2D)
  _compatible(into, from) ||
    throw(ArgumentError("merge!: 2D histograms must share box, cell counts, and lag"))
  into.counts .+= from.counts
  into.powersums .+= from.powersums
  return into
end

end # module Histograms
