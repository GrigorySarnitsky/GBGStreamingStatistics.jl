# Histograms — the carriers (universal lower half): the finely binned empirical measures that are
# the canonical summary `T`. Linear (WARPing) binning, per-bin counts, per-bin standardized central
# moments, and per-bin power-sums. Rationale: `docs/ARCHITECTURE.md` §3.

module Histograms

using ..Accumulators: AbstractAccumulator

export Histogram1D, Histogram2D, bin_index, observe_pairs!


"""
    Histogram1D <: AbstractAccumulator

A **finely binned empirical measure** of a conditioning variable ``X`` — the canonical summary
``T(\\hat\\mu_n)`` (`docs/ARCHITECTURE.md` §3). Linear (WARPing) bins over a fixed range (default
``\\pm 5\\sigma``, ``\\gtrsim 10^3`` bins), so the bin width ``b \\approx \\sigma/100 \\ll h`` makes
the summary *approximately sufficient*: the ``O(b^2)`` binning error is dominated by the irreducible
``O(h^2)`` smoothing bias.

Each bin carries its **count** and **per-bin standardized central moments to order ``\\approx 12``**
— the heavy-tailed channels need the high order (the Pawula ``Q``-criterion's ``M_6`` has a sampling
variance reaching the 12th moment) — accumulated by the Pébay + Kahan recursion. Per-bin
**power-sums** ``S_{rs} = \\sum_{i \\in \\text{bin}} Y_i^{r}\\,(X_i - c)^{s}`` (about the bin center
``c``) are the local sufficient statistic for local-polynomial regression. Counts are integer; moment
accumulators are `Float64` by default (accumulate wider than the data).

!!! note "Scaffold stage"
    The fields below are a provisional skeleton; storage is a generic `AbstractArray` (no scalar
    indexing of device arrays) so the carrier is GPU-resident. Not yet implemented.
"""
struct Histogram1D{T,C<:AbstractVector,A<:AbstractArray} <: AbstractAccumulator
  lo::T          # lower edge of the binned range
  hi::T          # upper edge of the binned range
  counts::C      # per-bin integer counts (UInt64), a generic array for GPU residency
  moments::A     # per-bin standardized central moments to order ≈ 12 (Float64)
  powersums::A   # per-bin power-sums S_rs, the local sufficient statistic for regression
end


"""
    Histogram2D <: AbstractAccumulator

A **joint binned empirical measure** of ``(X, Y)`` — one per lag ``\\tau`` — for conditional and
two-time transition densities (`docs/ARCHITECTURE.md` §3, §4.3), with the same per-bin carriers as
[`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) on the product grid. For a two-time quantity the conditioning value ``X(t)``
and the response ``Y(t+\\tau)`` must be paired *per source*: the entry is made at the response step
into the bin of ``X(t)`` (two separate marginals cannot recover the joint — for the
Ornstein–Uhlenbeck oracle both marginals are the same stationary Gaussian and the entire
``\\tau``-dependence lives in the pair). The lag ``\\tau`` is an explicit label on the carrier.

!!! note "Scaffold stage"
    Provisional skeleton; not yet implemented.
"""
struct Histogram2D{T,C<:AbstractArray,A<:AbstractArray} <: AbstractAccumulator
  xlo::T         # lower edge of the X range
  xhi::T         # upper edge of the X range
  ylo::T         # lower edge of the Y range
  yhi::T         # upper edge of the Y range
  lag::Int       # the lag τ this joint carrier is accumulated at (0 for the static joint)
  counts::C      # per-cell integer counts on the (X, Y) product grid
  moments::A     # per-cell standardized central moments (Float64)
end


"""
    bin_index(summary::AbstractAccumulator, x) -> i

The **binning map**: the index ``i`` of the linear (WARPing) bin containing the conditioning value
``x``. Exposed so a consumer can pre-quantize a snapshot buffer (e.g. store compact `UInt16` bin
indices when memory is tight) before forming lagged pairs for [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!).

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function bin_index end


"""
    observe_pairs!(summary, x_past, y_now; lag) -> summary

Fold a batch of **pre-formed, lag-labelled pairs** ``\\bigl(X(t),\\, Y(t+\\tau)\\bigr)`` into the
per-lag carrier in place — the memoryless surface for two-time statistics (conditional and transition
densities, conditional increment moments; `docs/ARCHITECTURE.md` §3, §4.3–4.4). Everything that
*produces* the pairs — the snapshot buffer of past ``X(t)``, which sources to sample, which lags to
form, the origin cadence — is the **consumer's** responsibility; this package only bins the pairs it
is handed, keeping the small per-lag histograms and never the ``O(ML)`` snapshot buffer.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function observe_pairs! end

end # module Histograms
