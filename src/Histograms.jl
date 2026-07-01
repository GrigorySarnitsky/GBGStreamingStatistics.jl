# Histograms — the summaries (universal lower half): the finely binned empirical measures that are the
# canonical summary `T`. Linear (WARPing) binning; per-bin counts + power-sums about the bin center,
# accumulated by one portable KernelAbstractions + Atomix scatter kernel. Rationale: `docs/ARCHITECTURE.md` §3.

module Histograms

import ..Accumulators: observe!, scale!
import Base: merge!
using ..Accumulators: AbstractAccumulator
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
import Atomix

export Histogram1D, Histogram2D, bin_index, observe_pairs!, scale!


"""
    Histogram1D <: AbstractAccumulator
    Histogram1D(lower, upper, nbins; order = 4, accumulator = Float64, count = UInt64)

A **finely binned empirical measure** of a conditioning variable ``X`` — the canonical summary
``T(\\hat\\mu_n)``. `nbins` equal-width (WARPing) bins tile the half-open range ``[\\texttt{lower}, \\texttt{upper})``;
choosing the range to cover the data (e.g. ``\\pm 5\\sigma`` for the data's standard deviation
``\\sigma``) with ``\\gtrsim 10^3`` bins makes the bin width ``b \\approx \\sigma/100`` far below any
post-hoc smoothing bandwidth ``h`` (``b \\ll h``), so the summary is *approximately sufficient*.

Each bin carries its **count** and the **power-sums** ``S_p = \\sum_{i \\in \\text{bin}}(X_i - c_j)^p``
about the fixed bin center ``c_j``, for ``p = 1 \\dots \\texttt{order}``. Because ``X_i - c_j`` is small within
a bin these are well-conditioned, so accumulation is a plain atomic scatter-add (`observe!`) that runs
identically on CPU and GPU, and `merge!` is elementwise addition. The power-sums use the *actual*
``X_i`` (not the center), so the global central moments recovered post hoc (by the shift formula) are
exact, and stable for the even orders (the odd orders, ``\\approx 0``, keep weaker relative precision).
The default `order = 4` recovers up to ``M_4``; an even higher moment such as ``M_6`` needs `order = 6`.
Counts are `UInt64` by default (the additive monoid) and power-sums are accumulated in `accumulator`
(`Float64` by default) regardless of the input element type — *accumulate wider than the data*. Set
`count = Float64` for the **Float-count variant** that [`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!)
and weighted [`merge!`](@ref Base.merge!) require — the cone the hyperstatistical parameter-mixture and
exponential forgetting are built from (a fractional weighted count is not an integer).

Samples outside ``[\\texttt{lower}, \\texttt{upper})`` are **skipped** (not counted) — set the range to cover the data.

# Fields

The summary's whole state, each reachable by dot — `summary.lower`, `summary.scale`,
`summary.counts`, … — for serialization, inspection, or building one by hand:
- `lower`, `upper` — the range edges (`upper > lower`).
- `scale` — the precomputed binning factor `nbins / (upper - lower)`, stored (not recomputed) so that
  `bin_index` is a single multiply.
- `counts` — per-bin counts, length `nbins` (`UInt64`, or the `count` type; `nbins == length(counts)`).
- `powersums` — per-bin power-sums, ``\\texttt{order} \\times \\texttt{nbins}`` (row ``p`` is ``S_p``; `order == size(powersums, 1)`).

A freshly constructed summary is **empty**: `counts` and `powersums` are all zeros (the `merge!`
identity), and you `observe!` data into it. `merge!` adds `counts` and `powersums` elementwise, so
summaries on the same binning combine by saving and adding these two arrays. The recovered moments are
read out with `mean`, `var`, `std`, and `central_moment`.
"""
struct Histogram1D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  lower::T          # lower edge of the binned range
  upper::T          # upper edge of the binned range
  scale::T       # nbins / (upper − lower); bin of x is floor((x − lower)·scale) + 1
  counts::C      # length nbins, per-bin integer counts (UInt64)
  powersums::M   # order × nbins, per-bin power-sums S_1 … S_order about the bin center
end


function Histogram1D(lower::Real, upper::Real, nbins::Integer;
                     order::Integer = 4, accumulator::Type{<:AbstractFloat} = Float64,
                     count::Type{<:Real} = UInt64)
  upper > lower || throw(ArgumentError("Histogram1D requires upper > lower (got lower = $lower, upper = $upper)"))
  nbins ≥ 1 || throw(ArgumentError("Histogram1D requires nbins ≥ 1 (got $nbins)"))
  order ≥ 2 || throw(ArgumentError("Histogram1D requires order ≥ 2 (got $order)"))
  lowerT, upperT = promote(float(lower), float(upper))
  scale = oftype(lowerT, nbins) / (upperT - lowerT)
  counts = zeros(count, nbins)
  powersums = zeros(accumulator, order, nbins)
  return Histogram1D(lowerT, upperT, scale, counts, powersums)
end


_nbins(h::Histogram1D) = length(h.counts)
_order(h::Histogram1D) = size(h.powersums, 1)
_bincenter(h::Histogram1D, i::Integer) = h.lower + (oftype(h.lower, i) - oftype(h.lower, 0.5)) / h.scale
_total(h::Histogram1D) = sum(h.counts)


"""
    bin_index(summary, x) -> i

The **binning map**: the index ``i`` of the bin containing `x`. For a
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) over ``[\\texttt{lower}, \\texttt{upper})`` with `nbins`
bins, ``i = \\lfloor (x - \\texttt{lower})\\,\\texttt{scale} \\rfloor + 1``; a valid bin satisfies
``1 \\le i \\le \\texttt{nbins}``, and `x` outside the range yields ``i < 1`` or ``i > \\texttt{nbins}``. Exposed so a
consumer can pre-quantize a snapshot buffer (e.g. store compact `UInt16` bin indices) before forming
lagged pairs for [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!).
"""
@inline bin_index(h::Histogram1D, x) = floor(Int, (x - h.lower) * h.scale) + 1


# The portable scatter-add: one work-item per sample atomically accumulates the count and the
# per-bin power-sums about the bin center, directly into the preallocated summary. Runs on CPU and
# any GPU backend (KernelAbstractions + Atomix). Out-of-range samples are skipped.
@kernel function _scatter1d!(counts, powersums, @Const(samples), lower, scale, nbins, order)
  s = @index(Global, Linear)
  if s ≤ length(samples)
    @inbounds begin
      x = samples[s]
      i = floor(Int, (x - lower) * scale) + 1
      if 1 ≤ i ≤ nbins
        Atomix.@atomic counts[i] += one(eltype(counts))
        center = lower + (oftype(lower, i) - oftype(lower, 0.5)) / scale
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
  _scatter1d!(backend)(h.counts, h.powersums, samples, h.lower, h.scale, _nbins(h), _order(h);
                       ndrange = length(samples))
  synchronize(backend)
  return h
end


# Scalar fold — a direct host update (scalar-indexes the summary, so host arrays only; batch on GPU).
function observe!(h::Histogram1D, x::Real)
  i = floor(Int, (x - h.lower) * h.scale) + 1
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
  a.lower == b.lower && a.upper == b.upper && _nbins(a) == _nbins(b) && _order(a) == _order(b)


# scale! and a weighted merge! produce fractional counts, so they need a Float-count summary
# (`count = Float64`); the default UInt64 counts would truncate a weighted count w·n. Guard, clearly.
_require_float_counts(h, op) =
  eltype(h.counts) <: AbstractFloat || throw(ArgumentError(
    "$op needs a Float-count summary: construct with `count = Float64` (the default `UInt64` counts cannot hold a fractional weighted count)"))


function merge!(into::Histogram1D, from::Histogram1D; weight = 1)
  _compatible(into, from) ||
    throw(ArgumentError("merge!: histograms must share lower, upper, bin count, and order"))
  if isone(weight)
    into.counts .+= from.counts          # counts and power-sums add (same binning ⇒ same centers)
    into.powersums .+= from.powersums
  else
    _require_float_counts(into, "a weighted merge!")
    into.counts .+= weight .* from.counts          # the cone: into += weight·from
    into.powersums .+= weight .* from.powersums
  end
  return into
end


# The cone scalar action: reweight every count and power-sum by `weight` in place (hyperstatistics /
# exponential forgetting). Float counts only — an integer count cannot hold a fractional weight·n.
function scale!(h::Histogram1D, weight)
  _require_float_counts(h, "scale!")
  h.counts .*= weight
  h.powersums .*= weight
  return h
end


"""
    Histogram2D <: AbstractAccumulator
    Histogram2D(xlower, xupper, nx, ylower, yupper, ny; lag = 0,
                accumulator = Float64, count = UInt64)

A **joint binned empirical measure** of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities. ``\\texttt{nx} \\times \\texttt{ny}`` equal-width cells tile
``[\\texttt{xlower}, \\texttt{xupper}) \\times [\\texttt{ylower}, \\texttt{yupper})``; each cell `(ix, iy)` (flattened `ix + (iy-1)*nx`) carries its
**count** and the order-2 bivariate **power-sums** ``S_{rs} = \\sum_{(X,Y)\\in\\text{cell}} (X -
c^x)^r (Y - c^y)^s`` about the cell center ``(c^x, c^y)`` — ``S_{10}, S_{01}, S_{20}, S_{11},
S_{02}`` — the local sufficient statistic for the joint density and for local-linear conditioning.
For a two-time quantity ``X(t)`` and the response ``Y(t+\\tau)`` are paired *per source*; the lag
``\\tau`` is an explicit label on the summary. Pairs outside the box are skipped. Same atomic-scatter
accumulation, `merge!`-by-addition, `Float64` accumulation, and `count = Float64` Float-count / cone
variant (for [`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!) and weighted
[`merge!`](@ref Base.merge!)) as [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D).

# Fields

The summary's whole state, each reachable by dot (`summary.xlower`, `summary.xscale`, `summary.lag`, …):
- `xlower`, `xupper`, `ylower`, `yupper`, `nx`, `ny` — the box and cell counts.
- `xscale`, `yscale` — the precomputed binning factors `nx/(xupper-xlower)`, `ny/(yupper-ylower)`, stored.
- `lag` — the lag ``\\tau`` this joint summary is accumulated at (`0` for the static joint).
- `counts` — per-cell counts, length `nx*ny` (cell `(ix, iy)` at `ix + (iy-1)*nx`; `UInt64`, or the `count` type).
- `powersums` — ``5 \\times (\\texttt{nx} \\cdot \\texttt{ny})``: rows ``S_{10}, S_{01}, S_{20}, S_{11}, S_{02}`` about the cell center.

A freshly constructed summary is **empty** (`counts`, `powersums` all zeros). `merge!` adds `counts`
and `powersums` (requires matching box, cell counts, and `lag`).
"""
struct Histogram2D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  xlower::T; xupper::T; xscale::T     # X range and nx/(xupper−xlower)
  ylower::T; yupper::T; yscale::T     # Y range and ny/(yupper−ylower)
  nx::Int
  ny::Int
  lag::Int                       # the lag τ this joint summary is accumulated at (0 for the static joint)
  counts::C                      # length nx·ny, per-cell integer counts (UInt64)
  powersums::M                   # 5 × (nx·ny): rows S_10, S_01, S_20, S_11, S_02 about the cell center
end


function Histogram2D(xlower::Real, xupper::Real, nx::Integer, ylower::Real, yupper::Real, ny::Integer;
                     lag::Integer = 0, accumulator::Type{<:AbstractFloat} = Float64,
                     count::Type{<:Real} = UInt64)
  (xupper > xlower && yupper > ylower) ||
    throw(ArgumentError("Histogram2D requires xupper > xlower and yupper > ylower"))
  (nx ≥ 1 && ny ≥ 1) || throw(ArgumentError("Histogram2D requires nx ≥ 1 and ny ≥ 1"))
  xlowerT, xupperT, ylowerT, yupperT = promote(float(xlower), float(xupper), float(ylower), float(yupper))
  T = typeof(xlowerT)
  xscale = T(nx) / (xupperT - xlowerT)
  yscale = T(ny) / (yupperT - ylowerT)
  counts = zeros(count, nx * ny)
  powersums = zeros(accumulator, 5, nx * ny)
  return Histogram2D(xlowerT, xupperT, xscale, ylowerT, yupperT, yscale, Int(nx), Int(ny), Int(lag), counts, powersums)
end


_ncells(h::Histogram2D) = length(h.counts)
_total(h::Histogram2D) = sum(h.counts)
_xbincenter(h::Histogram2D, ix::Integer) = h.xlower + (oftype(h.xlower, ix) - oftype(h.xlower, 0.5)) / h.xscale
_ybincenter(h::Histogram2D, iy::Integer) = h.ylower + (oftype(h.ylower, iy) - oftype(h.ylower, 0.5)) / h.yscale


# Per-X-column reduction of a `Histogram2D` (the regression-ready sufficient statistic): for each X-bin,
# the count ``N``, the response total ``\sum Y`` and ``\sum Y^2`` — recovered *exactly* from the per-cell
# counts and power-sums (``S_{01}`` / ``S_{02}`` use the actual ``Y``, not the cell center). Host-side;
# a no-op `Array` for CPU summaries. Shared by the conditional-mean / increment-moment / bandwidth consumers.
function _column_sums(summary::Histogram2D)
  counts = Array(summary.counts)
  S = Array(summary.powersums)
  nx = summary.nx
  ny = summary.ny
  cx = [_xbincenter(summary, ix) for ix in 1:nx]
  N = zeros(Float64, nx)
  sumY = zeros(Float64, nx)
  sumY2 = zeros(Float64, nx)
  for iy in 1:ny
    cy = _ybincenter(summary, iy)
    for ix in 1:nx
      cell = ix + (iy - 1) * nx
      nc = counts[cell]
      iszero(nc) && continue
      s01 = S[2, cell]
      s02 = S[5, cell]
      N[ix] += nc
      sumY[ix] += cy * nc + s01                       # Σ Y over the column, exact
      sumY2[ix] += cy * cy * nc + 2 * cy * s01 + s02  # Σ Y² over the column, exact
    end
  end
  return cx, N, sumY, sumY2
end


# The 2D scatter-add: one work-item per pair atomically accumulates the per-cell count and the five
# order-2 bivariate power-sums about the cell center. Out-of-box pairs are skipped.
@kernel function _scatter2d!(counts, powersums, @Const(xs), @Const(ys),
                             xlower, xscale, nx, ylower, yscale, ny)
  s = @index(Global, Linear)
  if s ≤ length(xs)
    @inbounds begin
      x = xs[s]
      y = ys[s]
      ix = floor(Int, (x - xlower) * xscale) + 1
      iy = floor(Int, (y - ylower) * yscale) + 1
      if 1 ≤ ix ≤ nx && 1 ≤ iy ≤ ny
        cell = ix + (iy - 1) * nx
        Atomix.@atomic counts[cell] += one(eltype(counts))
        dx = x - (xlower + (oftype(xlower, ix) - oftype(xlower, 0.5)) / xscale)
        dy = y - (ylower + (oftype(ylower, iy) - oftype(ylower, 0.5)) / yscale)
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
place (`x_past` holds ``X(t)``, `y_now` holds ``Y(t+\\tau)``, equal length) — the memoryless surface for
two-time statistics. The pairs are binned by an atomic
scatter-add (one work-item per pair, CPU and GPU alike); everything that *produces* the pairs — the
snapshot buffer of past ``X(t)``, which sources to sample, which lags to form, the origin cadence — is
the **consumer's** responsibility, formed at the summary's lag ``\\tau``.
"""
function observe_pairs!(h::Histogram2D, x_past, y_now)
  length(x_past) == length(y_now) ||
    throw(ArgumentError("observe_pairs!: x_past and y_now must have equal length"))
  backend = get_backend(h.counts)
  _scatter2d!(backend)(h.counts, h.powersums, x_past, y_now,
                       h.xlower, h.xscale, h.nx, h.ylower, h.yscale, h.ny; ndrange = length(x_past))
  synchronize(backend)
  return h
end


_compatible(a::Histogram2D, b::Histogram2D) =
  a.xlower == b.xlower && a.xupper == b.xupper && a.ylower == b.ylower && a.yupper == b.yupper &&
  a.nx == b.nx && a.ny == b.ny && a.lag == b.lag

function merge!(into::Histogram2D, from::Histogram2D; weight = 1)
  _compatible(into, from) ||
    throw(ArgumentError("merge!: 2D histograms must share box, cell counts, and lag"))
  if isone(weight)
    into.counts .+= from.counts
    into.powersums .+= from.powersums
  else
    _require_float_counts(into, "a weighted merge!")
    into.counts .+= weight .* from.counts          # the cone: into += weight·from
    into.powersums .+= weight .* from.powersums
  end
  return into
end


function scale!(h::Histogram2D, weight)
  _require_float_counts(h, "scale!")
  h.counts .*= weight
  h.powersums .*= weight
  return h
end

end # module Histograms
