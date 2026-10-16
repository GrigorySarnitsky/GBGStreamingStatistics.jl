# Histograms — the summaries (universal lower half): the finely binned empirical measures that are the
# canonical summary `T`. Linear (WARPing) binning; per-bin counts + power-sums about the bin center,
# accumulated by one portable KernelAbstractions + Atomix scatter kernel. Rationale: `docs/ARCHITECTURE.md` §3.

module Histograms

import ..Accumulators: observe!, scale!
import Base: merge!
using ..Accumulators: AbstractAccumulator
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
import Atomix

export Histogram1D, Histogram2D, observe_pairs!, scale!, write_histogram, read_histogram


"""
    Histogram1D <: AbstractAccumulator
    Histogram1D(lower, upper, nbins; order = 4, accumulator_type = Float64,
                count_type = UInt64)

A **finely binned empirical measure** of a scalar variate ``X`` — the canonical summary
``T(\\hat\\mu_n)``. The `nbins` equal-width bins (WARPing — weighted averaging of rounded points) tile the half-open range ``[\\texttt{lower}, \\texttt{upper})``;
choosing the range to cover the data (e.g. ``\\pm 5\\sigma`` for the data's standard deviation
``\\sigma``) with ``\\gtrsim 10^3`` bins makes the bin width ``b \\approx \\sigma/100`` far below any
post-hoc smoothing bandwidth ``h`` (``b \\ll h``), so the summary is *approximately sufficient*.

Each bin carries its **count** and the **power-sums** ``S_p = \\sum_{i \\in \\text{bin}}(X_i - c_j)^p``
about the fixed bin center ``c_j``, for ``p = 1 \\dots \\texttt{order}``. Because ``X_i - c_j`` is small within
a bin these are well-conditioned, so accumulation is a plain atomic scatter-add
([`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!)) that runs
identically on CPU and GPU, and `merge!` is elementwise addition. The power-sums use the *actual*
``X_i`` (not the center), so the global central moments recovered post hoc (by the shift formula — the binomial re-centering
from the bin center to the sample mean) are exact, and stable for the even orders (the odd orders, ``\\approx 0``, keep weaker relative precision).
The default `order = 4` recovers up to ``M_4``; an even higher moment such as ``M_6`` needs `order = 6`.

Counts are `UInt64` by default (summaries add) and power-sums are accumulated in `accumulator_type`
(`Float64` by default) regardless of the input element type — *accumulate wider than the data*. Set
`count_type = Float64` for the **Float-count variant** that [`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!)
and weighted [`merge!`](@ref Base.merge!) require — the convex cone (summaries add and scale) that
post-hoc re-weighting of a family of summaries (hyperstatistics) and down-weighting old data
(exponential forgetting) need (a fractional weighted count is not an integer).

Samples outside ``[\\texttt{lower}, \\texttt{upper})`` are **skipped** (not counted).

# Arguments

- `lower`, `upper` — the edges of the half-open range ``[\\texttt{lower}, \\texttt{upper})`` (`upper > lower`).
- `nbins` — the number of equal-width bins.
- `order` — the highest power-sum order ``p`` stored per bin (default `4`, recovering up to ``M_4``).
- `accumulator_type` — the element type the power-sums are accumulated in (default `Float64`).
- `count_type` — the element type of the counts: `UInt64` (the default) or `Float64` for the Float-count
  variant.


A freshly constructed summary is **empty**: `counts` and `powersums` are all zeros (the `merge!`
identity), and you `observe!` data into it; `merge!` adds `counts` and `powersums` elementwise, so
same-binning summaries combine by adding these two arrays. The recovered moments are
read out with [`central_moment`](@ref GBGStreamingStatistics.Expectations.central_moment) and, for the
first two, the `Statistics` extensions `mean`, `var`, and `std`.

```jldoctest
julia> summary = Histogram1D(-5.0, 5.0, 1000);       # 1000 bins over [-5, 5)

julia> observe!(summary, [0.3, -1.2, 0.3, 2.1]);     # fold in a stream, never stored

julia> Int(sum(summary.counts))                      # the four in-range samples counted
4
```

# Fields

The fields hold the summary's whole **internal state**, distinct from the constructor arguments (`scale`
is derived, `counts` and `powersums` allocated); each is reachable by dot (`summary.lower`,
`summary.scale`, `summary.counts`, …), for serialization, inspection, or building one by hand:
- `lower`, `upper` — the range edges (`upper > lower`).
- `scale` — the precomputed binning factor `nbins / (upper - lower)`, stored (not recomputed) so that
  the binning map is a single multiply.
- `counts` — per-bin counts, length `nbins` (`UInt64`, or `count_type`; `nbins == length(counts)`).
- `powersums` — per-bin power-sums, ``\\texttt{order} \\times \\texttt{nbins}`` (row ``p`` is ``S_p``; `order == size(powersums, 1)`).
"""
struct Histogram1D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  lower::T          # lower edge of the binned range
  upper::T          # upper edge of the binned range
  scale::T       # nbins / (upper − lower); bin of x is floor((x − lower)·scale) + 1
  counts::C      # length nbins, per-bin integer counts (UInt64)
  powersums::M   # order × nbins, per-bin power-sums S_1 … S_order about the bin center
end


function Histogram1D(lower::Real, upper::Real, nbins::Integer;
                     order::Integer = 4, accumulator_type::Type{<:AbstractFloat} = Float64,
                     count_type::Type{<:Real} = UInt64)
  upper > lower || throw(ArgumentError("Histogram1D requires upper > lower (got lower = $lower, upper = $upper)"))
  nbins ≥ 1 || throw(ArgumentError("Histogram1D requires nbins ≥ 1 (got $nbins)"))
  order ≥ 2 || throw(ArgumentError("Histogram1D requires order ≥ 2 (got $order)"))
  lowerT, upperT = promote(float(lower), float(upper))
  scale = oftype(lowerT, nbins) / (upperT - lowerT)
  counts = zeros(count_type, nbins)
  powersums = zeros(accumulator_type, order, nbins)
  return Histogram1D(lowerT, upperT, scale, counts, powersums)
end


_nbins(h::Histogram1D) = length(h.counts)
_order(h::Histogram1D) = size(h.powersums, 1)
_bincenter(h::Histogram1D, i::Integer) = h.lower + (oftype(h.lower, i) - oftype(h.lower, 0.5)) / h.scale
_total(h::Histogram1D) = sum(h.counts)


# `bin_index` (the binning map, `floor(Int, (x - lower) * scale) + 1`) is **archived**: it had no internal
# caller (the scatter kernels inline the map), and the public path it was exposed for — pre-quantizing the
# snapshot buffer to compact bin indices and folding them through an index-taking `observe_pairs!` (lossless,
# since the conditioning readers use only the x-cell, never its sub-bin offset) — is deferred. Restore and
# re-export it with that path.
# @inline bin_index(h::Histogram1D, x) = floor(Int, (x - h.lower) * h.scale) + 1


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


function observe!(summary::Histogram1D, batch)
  backend = get_backend(summary.counts)
  _scatter1d!(backend)(summary.counts, summary.powersums, batch, summary.lower, summary.scale,
                       _nbins(summary), _order(summary); ndrange = length(batch))
  synchronize(backend)
  return summary
end


# Scalar fold — a direct host update (scalar-indexes the summary, so host arrays only; batch on GPU).
function observe!(summary::Histogram1D, batch::Real)
  i = floor(Int, (batch - summary.lower) * summary.scale) + 1
  (1 ≤ i ≤ _nbins(summary)) || return summary
  @inbounds summary.counts[i] += one(eltype(summary.counts))
  dx = batch - _bincenter(summary, i)
  dxp = dx
  @inbounds for p in 1:_order(summary)
    summary.powersums[p, i] += dxp
    dxp *= dx
  end
  return summary
end


_compatible(a::Histogram1D, b::Histogram1D) =
  a.lower == b.lower && a.upper == b.upper && _nbins(a) == _nbins(b) && _order(a) == _order(b)


# scale! and a weighted merge! produce fractional counts, so they need a Float-count summary
# (`count_type = Float64`); the default UInt64 counts would truncate a weighted count w·n. Guard, clearly.
_require_float_counts(h, op) =
  eltype(h.counts) <: AbstractFloat || throw(ArgumentError(
    "$op needs a Float-count summary: construct with `count_type = Float64` (the default `UInt64` counts cannot hold a fractional weighted count)"))


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
function scale!(summary::Histogram1D, weight)
  _require_float_counts(summary, "scale!")
  summary.counts .*= weight
  summary.powersums .*= weight
  return summary
end


"""
    Histogram2D <: AbstractAccumulator
    Histogram2D(xlower, xupper, nx, ylower, yupper, ny; lag = 0,
                accumulator_type = Float64, count_type = UInt64)

A **joint binned empirical measure** of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities. Its ``\\texttt{nx} \\times \\texttt{ny}`` equal-width cells tile
``[\\texttt{xlower}, \\texttt{xupper}) \\times [\\texttt{ylower}, \\texttt{yupper})``; each cell `(ix, iy)` (flattened `ix + (iy-1)*nx`) carries its
**count** and the order-2 bivariate **power-sums** ``S_{rs} = \\sum_{(X,Y)\\in\\text{cell}} (X -
c^x)^r (Y - c^y)^s`` about the cell center ``(c^x, c^y)`` — ``S_{10}, S_{01}, S_{20}, S_{11},
S_{02}`` — the local sufficient statistic for the joint density and for local-linear conditioning.

For a two-time quantity, the value ``X(t)`` and the later response ``Y(t+\\tau)`` of one source (a stream
trajectory) are paired; the lag ``\\tau`` is an explicit label on the summary. Pairs outside
the box are skipped. Accumulation (atomic scatter), `merge!` by addition, `Float64` accumulation and the
`count_type = Float64` Float-count / cone variant (for
[`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!) and weighted [`merge!`](@ref Base.merge!))
are as for [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D).

# Arguments

- `xlower`, `xupper`, `ylower`, `yupper` — the edges of the box
  ``[\\texttt{xlower}, \\texttt{xupper}) \\times [\\texttt{ylower}, \\texttt{yupper})``.
- `nx`, `ny` — the numbers of equal-width cells in ``x`` and in ``y``.
- `lag` — the lag ``\\tau`` the pairs are formed at (default `0`, the static joint).
- `accumulator_type` — the element type the power-sums are accumulated in (default `Float64`).
- `count_type` — the element type of the counts: `UInt64` (the default) or `Float64` for the Float-count
  variant.


A freshly constructed summary is **empty** (`counts`, `powersums` all zeros), and `merge!` adds `counts`
and `powersums` (requires matching box, cell counts, and `lag`).

```jldoctest
julia> joint = Histogram2D(-5.0, 5.0, 100, -5.0, 5.0, 100);   # 100×100 cells over the box

julia> observe_pairs!(joint, [0.3, -1.2], [0.5, 0.4]);        # (X, Y) pairs, paired per source

julia> Int(sum(joint.counts))                                 # the two in-box pairs counted
2
```

# Fields

The fields hold the summary's whole **internal state**, distinct from the constructor arguments
(`xscale`/`yscale` derived, `counts` and `powersums` allocated); each is reachable by dot
(`summary.xlower`, `summary.xscale`, `summary.lag`, …):
- `xlower`, `xupper`, `ylower`, `yupper`, `nx`, `ny` — the box and cell counts.
- `xscale`, `yscale` — the precomputed binning factors `nx/(xupper-xlower)`, `ny/(yupper-ylower)`, stored.
- `lag` — the lag ``\\tau`` this joint summary is accumulated at (`0` for the static joint).
- `counts` — per-cell counts, length `nx*ny` (cell `(ix, iy)` at `ix + (iy-1)*nx`; `UInt64`, or `count_type`).
- `powersums` — ``5 \\times (\\texttt{nx} \\cdot \\texttt{ny})``: rows ``S_{10}, S_{01}, S_{20}, S_{11}, S_{02}`` about the cell center.
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
                     lag::Integer = 0, accumulator_type::Type{<:AbstractFloat} = Float64,
                     count_type::Type{<:Real} = UInt64)
  (xupper > xlower && yupper > ylower) ||
    throw(ArgumentError("Histogram2D requires xupper > xlower and yupper > ylower"))
  (nx ≥ 1 && ny ≥ 1) || throw(ArgumentError("Histogram2D requires nx ≥ 1 and ny ≥ 1"))
  xlowerT, xupperT, ylowerT, yupperT = promote(float(xlower), float(xupper), float(ylower), float(yupper))
  T = typeof(xlowerT)
  xscale = T(nx) / (xupperT - xlowerT)
  yscale = T(ny) / (yupperT - ylowerT)
  counts = zeros(count_type, nx * ny)
  powersums = zeros(accumulator_type, 5, nx * ny)
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

Fold a batch of **pre-formed pairs** ``\\bigl(X(t),  Y(t+\\tau)\\bigr)`` into the per-lag summary in
place (`x_past` holds ``X(t)``, `y_now` holds ``Y(t+\\tau)``, equal length) — the entry point that stores
no sample. The pairs are binned by an atomic scatter-add (one work-item per pair, CPU or GPU);
everything that *produces* the pairs — the buffer of past ``X(t)`` snapshots, which sources (stream
trajectories) to sample, which lags to form, how often to start a pair — is the **consumer's** (the
calling code's) responsibility, formed at the summary's lag ``\\tau``.

# Arguments

- `summary` — the per-lag joint `Histogram2D` the pairs are folded into.
- `x_past` — the conditioning values ``X(t)``, one per pair.
- `y_now` — the response values ``Y(t+\\tau)``, one per pair, in the same order as `x_past`.
"""
function observe_pairs!(summary::Histogram2D, x_past, y_now)
  length(x_past) == length(y_now) ||
    throw(ArgumentError("observe_pairs!: x_past and y_now must have equal length"))
  backend = get_backend(summary.counts)
  _scatter2d!(backend)(summary.counts, summary.powersums, x_past, y_now,
                       summary.xlower, summary.xscale, summary.nx, summary.ylower, summary.yscale, summary.ny; ndrange = length(x_past))
  synchronize(backend)
  return summary
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


function scale!(summary::Histogram2D, weight)
  _require_float_counts(summary, "scale!")
  summary.counts .*= weight
  summary.powersums .*= weight
  return summary
end


# --- HDF5 persistence: bare generics declared here so `?write_histogram` and the docs work without HDF5
# present; the methods (on `HDF5.File` / `Group`) live in `ext/GBGStreamingStatisticsHDF5Ext.jl`, loaded
# automatically when the HDF5 package is. ---

"""
    write_histogram(parent, name, summary)

Persist a `Histogram1D` or `Histogram2D` `summary` as a **new HDF5 group** `name` under `parent` (an open
HDF5 `File` or `Group`), returning that group. The summary is the durable artifact of the streaming design
— a fixed-size object — so this is how you checkpoint it, resume from it, or combine summaries stored
across runs; the reload is exact.

The group is **self-describing**: the arrays become datasets (`counts`, `powersums`) and the scalar
parameters become attributes (`type`, `version`, the range edges, `nbins` / `nx` / `ny`, `order`, `lag`),
so any HDF5 tool can read it. Device (GPU) arrays are copied to the host on write, and the count element
type is preserved, so the `count_type = Float64` cone variant round-trips. The sugar `parent[name] = summary`
is equivalent.

The handle `parent` is one you opened and own — this never opens or closes a file. A package extension
provides this method; it is available only when `HDF5` is loaded.

# Arguments

- `parent` — an open HDF5 `File` or `Group` to create the new subgroup under.
- `name` — the name of the group to create.
- `summary` — the `Histogram1D` or `Histogram2D` to write.
"""
function write_histogram end


"""
    read_histogram(parent, name)

Reconstruct the summary written by
[`write_histogram`](@ref GBGStreamingStatistics.Histograms.write_histogram) from the HDF5 group `name`
under `parent` (an open HDF5 `File` or `Group`), dispatching on the stored `type` attribute to return a
`Histogram1D` or a `Histogram2D`. The reload is **exact** — counts, power-sums, range, element types, and
(for a `Histogram2D`) the lag are recovered bit-for-bit, and the result `merge!`s with other summaries of
the same shape.

The reconstructed summary is a host (CPU) object regardless of where the written one lived; move it to a
device afterward if needed. A package extension provides this method; it is available only when `HDF5`
is loaded.

The arguments are `parent`, an open HDF5 `File` or `Group` holding the group, and `name`, the name (or
path) of the histogram group under `parent`.
"""
function read_histogram end

end # module Histograms
