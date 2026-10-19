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

A finely binned empirical measure of a scalar variate ``X`` — the canonical summary
``T(n\\mathbb{P}_n)``. The `nbins` equal-width bins (WARPing — weighted averaging of rounded points) tile the half-open range ``[\\texttt{lower}, \\texttt{upper})``;
choosing the range to cover the data (e.g. ``\\pm 5\\sigma_X`` for the data's standard deviation
``\\sigma_X``) with ``\\gtrsim 10^3`` bins makes the bin width ``\\Delta \\approx \\sigma_X/100`` far below any
post-hoc smoothing bandwidth ``h`` (``\\Delta \\ll h``), so the summary is *approximately sufficient*.

Each bin ``j`` carries its count and the power-sums ``S_p^{(j)} = \\sum_{i \\in \\text{bin } j}(X_i - c_j)^p``
about the fixed bin center ``c_j``, for ``p = 1 \\dots \\texttt{order}``. Because ``X_i - c_j`` is small within
a bin these are well-conditioned, so accumulation is a plain atomic scatter-add
([`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!)) that runs
identically on CPU and GPU, and `merge!` is elementwise addition. The power-sums use the *actual*
``X_i`` (not the center), so the global central moments recovered post hoc (by the shift formula — the binomial re-centering
from the bin center to the sample mean) are exact, and stable for the even orders (the odd orders, ``\\approx 0``, keep weaker relative precision).
The default `order = 4` recovers the central moments up to ``\\hat\\mu_4``; an even higher moment such as ``\\hat\\mu_6`` needs `order = 6`.

Counts are `UInt64` by default (summaries add) and power-sums are accumulated in `accumulator_type`
(`Float64` by default) regardless of the input element type — *accumulate wider than the data*. Set
`count_type = Float64` for the Float-count variant that [`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!)
and weighted [`merge!`](@ref Base.merge!) require — the convex cone (summaries add and scale) that
post-hoc re-weighting of a family of summaries (hyperstatistics) and down-weighting old data
(exponential forgetting) need (a fractional weighted count is not an integer).

Samples outside ``[\\texttt{lower}, \\texttt{upper})`` are skipped (not counted).

# Arguments

- `lower`, `upper` — the edges of the half-open range ``[\\texttt{lower}, \\texttt{upper})`` (`upper > lower`).
- `nbins` — the number of equal-width bins.
- `order` — the highest power-sum order stored per bin (default `4`, recovering up to ``\\hat\\mu_4``).
- `accumulator_type` — the element type the power-sums are accumulated in (default `Float64`).
- `count_type` — the element type of the counts: `UInt64` (the default) or `Float64` for the Float-count
  variant.


A freshly constructed summary is empty: `counts` and `powersums` are all zeros (the `merge!`
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

The fields hold the summary's whole internal state, distinct from the constructor arguments (`scale`
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


# `bin_index` (the binning map, `floor(Int, (x - lower) * scale) + 1`) is archived: it had no internal
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
# (`count_type = Float64`); the default UInt64 counts would truncate a weighted count λ·n. Guard, clearly.
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
    Histogram2D(x_lower, x_upper, n_x, y_lower, y_upper, n_y; lag = 0,
                accumulator_type = Float64, count_type = UInt64)

A joint binned empirical measure of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities. Its ``\\texttt{n\\_x} \\times \\texttt{n\\_y}`` equal-width cells tile
``[\\texttt{x_lower}, \\texttt{x_upper}) \\times [\\texttt{y_lower}, \\texttt{y_upper})``; each cell `(i_x, i_y)` (flattened `i_x + (i_y-1)*n_x`) carries its
count and the order-2 bivariate power-sums ``S_{pp'} = \\sum_{(X,Y)\\in\\text{cell}} (X -
c^x)^p (Y - c^y)^{p'}`` about the cell center ``(c^x, c^y)`` — ``S_{10}, S_{01}, S_{20}, S_{11},
S_{02}`` — the local sufficient statistic for the joint density and for local-linear conditioning.

For a two-time quantity, the value ``X(t)`` and the later response ``Y(t+\\tau)`` of one source (a stream
trajectory) are paired; the lag ``\\tau`` is an explicit label on the summary. Pairs outside
the box are skipped. Accumulation (atomic scatter), `merge!` by addition, `Float64` accumulation and the
`count_type = Float64` Float-count / cone variant (for
[`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!) and weighted [`merge!`](@ref Base.merge!))
are as for [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D).

# Arguments

- `x_lower`, `x_upper`, `y_lower`, `y_upper` — the edges of the box
  ``[\\texttt{x_lower}, \\texttt{x_upper}) \\times [\\texttt{y_lower}, \\texttt{y_upper})``.
- `n_x`, `n_y` — the numbers of equal-width cells in ``x`` and in ``y``.
- `lag` — the lag ``\\tau`` the pairs are formed at (default `0`, the static joint).
- `accumulator_type` — the element type the power-sums are accumulated in (default `Float64`).
- `count_type` — the element type of the counts: `UInt64` (the default) or `Float64` for the Float-count
  variant.


A freshly constructed summary is empty (`counts`, `powersums` all zeros), and `merge!` adds `counts`
and `powersums` (requires matching box, cell counts, and `lag`).

```jldoctest
julia> joint = Histogram2D(-5.0, 5.0, 100, -5.0, 5.0, 100);   # 100×100 cells over the box

julia> observe_pairs!(joint, [0.3, -1.2], [0.5, 0.4]);        # (X, Y) pairs, paired per source

julia> Int(sum(joint.counts))                                 # the two in-box pairs counted
2
```

# Fields

The fields hold the summary's whole internal state, distinct from the constructor arguments
(`x_scale`/`y_scale` derived, `counts` and `powersums` allocated); each is reachable by dot
(`summary.x_lower`, `summary.x_scale`, `summary.lag`, …):
- `x_lower`, `x_upper`, `y_lower`, `y_upper`, `n_x`, `n_y` — the box and cell counts.
- `x_scale`, `y_scale` — the precomputed binning factors `n_x/(x_upper-x_lower)`, `n_y/(y_upper-y_lower)`, stored.
- `lag` — the lag ``\\tau`` this joint summary is accumulated at (`0` for the static joint).
- `counts` — per-cell counts, length `n_x*n_y` (cell `(i_x, i_y)` at `i_x + (i_y-1)*n_x`; `UInt64`, or `count_type`).
- `powersums` — ``5 \\times (\\texttt{n\\_x} \\cdot \\texttt{n\\_y})``: rows ``S_{10}, S_{01}, S_{20}, S_{11}, S_{02}`` about the cell center.
"""
struct Histogram2D{T<:AbstractFloat,C<:AbstractVector,M<:AbstractMatrix} <: AbstractAccumulator
  x_lower::T; x_upper::T; x_scale::T     # x range and n_x/(x_upper−x_lower)
  y_lower::T; y_upper::T; y_scale::T     # y range and n_y/(y_upper−y_lower)
  n_x::Int
  n_y::Int
  lag::Int                       # the lag τ this joint summary is accumulated at (0 for the static joint)
  counts::C                      # length n_x·n_y, per-cell integer counts (UInt64)
  powersums::M                   # 5 × (n_x·n_y): rows S_10, S_01, S_20, S_11, S_02 about the cell center
end


function Histogram2D(x_lower::Real, x_upper::Real, n_x::Integer, y_lower::Real, y_upper::Real, n_y::Integer;
                     lag::Integer = 0, accumulator_type::Type{<:AbstractFloat} = Float64,
                     count_type::Type{<:Real} = UInt64)
  (x_upper > x_lower && y_upper > y_lower) ||
    throw(ArgumentError("Histogram2D requires x_upper > x_lower and y_upper > y_lower"))
  (n_x ≥ 1 && n_y ≥ 1) || throw(ArgumentError("Histogram2D requires n_x ≥ 1 and n_y ≥ 1"))
  xlowerT, xupperT, ylowerT, yupperT = promote(float(x_lower), float(x_upper), float(y_lower), float(y_upper))
  T = typeof(xlowerT)
  x_scale = T(n_x) / (xupperT - xlowerT)
  y_scale = T(n_y) / (yupperT - ylowerT)
  counts = zeros(count_type, n_x * n_y)
  powersums = zeros(accumulator_type, 5, n_x * n_y)
  return Histogram2D(xlowerT, xupperT, x_scale, ylowerT, yupperT, y_scale, Int(n_x), Int(n_y), Int(lag), counts, powersums)
end


_ncells(h::Histogram2D) = length(h.counts)
_total(h::Histogram2D) = sum(h.counts)
_xbincenter(h::Histogram2D, i_x::Integer) = h.x_lower + (oftype(h.x_lower, i_x) - oftype(h.x_lower, 0.5)) / h.x_scale
_ybincenter(h::Histogram2D, i_y::Integer) = h.y_lower + (oftype(h.y_lower, i_y) - oftype(h.y_lower, 0.5)) / h.y_scale


# Per-x-column reduction of a `Histogram2D` (the regression-ready sufficient statistic): for each x-bin j,
# the count N_j, the response totals Σ Y and Σ Y² — recovered *exactly* from the per-cell
# counts and power-sums (S_01 / S_02 use the actual Y, not the cell center). Host-side;
# a no-op `Array` for CPU summaries. Shared by the conditional-mean / increment-moment / bandwidth consumers.
function _column_sums(summary::Histogram2D)
  counts = Array(summary.counts)
  S = Array(summary.powersums)
  n_x = summary.n_x
  n_y = summary.n_y
  c_x = [_xbincenter(summary, i_x) for i_x in 1:n_x]
  N = zeros(Float64, n_x)
  sum_Y = zeros(Float64, n_x)
  sum_Y2 = zeros(Float64, n_x)
  for i_y in 1:n_y
    c_y = _ybincenter(summary, i_y)
    for i_x in 1:n_x
      cell = i_x + (i_y - 1) * n_x
      nc = counts[cell]
      iszero(nc) && continue
      s01 = S[2, cell]
      s02 = S[5, cell]
      N[i_x] += nc
      sum_Y[i_x] += c_y * nc + s01                       # Σ Y over the column, exact
      sum_Y2[i_x] += c_y * c_y * nc + 2 * c_y * s01 + s02  # Σ Y² over the column, exact
    end
  end
  return c_x, N, sum_Y, sum_Y2
end


# The 2D scatter-add: one work-item per pair atomically accumulates the per-cell count and the five
# order-2 bivariate power-sums about the cell center. Out-of-box pairs are skipped.
@kernel function _scatter2d!(counts, powersums, @Const(xs), @Const(ys),
                             x_lower, x_scale, n_x, y_lower, y_scale, n_y)
  s = @index(Global, Linear)
  if s ≤ length(xs)
    @inbounds begin
      x = xs[s]
      y = ys[s]
      i_x = floor(Int, (x - x_lower) * x_scale) + 1
      i_y = floor(Int, (y - y_lower) * y_scale) + 1
      if 1 ≤ i_x ≤ n_x && 1 ≤ i_y ≤ n_y
        cell = i_x + (i_y - 1) * n_x
        Atomix.@atomic counts[cell] += one(eltype(counts))
        dx = x - (x_lower + (oftype(x_lower, i_x) - oftype(x_lower, 0.5)) / x_scale)
        dy = y - (y_lower + (oftype(y_lower, i_y) - oftype(y_lower, 0.5)) / y_scale)
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

Fold a batch of pre-formed pairs ``\\bigl(X(t),  Y(t+\\tau)\\bigr)`` into the per-lag summary in
place (`x_past` holds ``X(t)``, `y_now` holds ``Y(t+\\tau)``, equal length) — the entry point that stores
no sample. The pairs are binned by an atomic scatter-add (one work-item per pair, CPU or GPU);
everything that *produces* the pairs — the buffer of past ``X(t)`` snapshots, which sources (stream
trajectories) to sample, which lags to form, how often to start a pair — is the consumer's (the
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
                       summary.x_lower, summary.x_scale, summary.n_x, summary.y_lower, summary.y_scale, summary.n_y; ndrange = length(x_past))
  synchronize(backend)
  return summary
end


_compatible(a::Histogram2D, b::Histogram2D) =
  a.x_lower == b.x_lower && a.x_upper == b.x_upper && a.y_lower == b.y_lower && a.y_upper == b.y_upper &&
  a.n_x == b.n_x && a.n_y == b.n_y && a.lag == b.lag

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

Persist a `Histogram1D` or `Histogram2D` `summary` as a new HDF5 group `name` under `parent` (an open
HDF5 `File` or `Group`), returning that group. The summary is the durable artifact of the streaming design
— a fixed-size object — so this is how you checkpoint it, resume from it, or combine summaries stored
across runs; the reload is exact.

The group is self-describing: the arrays become datasets (`counts`, `powersums`) and the scalar
parameters become attributes (`type`, `version`, the range edges, `nbins` / `n_x` / `n_y`, `order`, `lag`),
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
`Histogram1D` or a `Histogram2D`. The reload is exact — counts, power-sums, range, element types, and
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
