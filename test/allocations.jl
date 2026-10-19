using Test
using GBGStreamingStatistics
using Random
import GBGCodingTools
using GBGCodingTools.TestHarness: @alloc_below
include("public_surface.jl")

# The hot path is allocation-free (invariant 4): on a preallocated histogram, observe!/merge! allocate
# no bytes (measured on the CPU, after a warm-up call to force compilation). This is the contract that
# keeps a long GPU run from wedging once the device kernel lands.
#
# The audit is DECLARED to GBGCodingTools (RULE.TEST.ALLOC_ZERO): `@alloc_below` warms up, measures, and
# records the operation as audited, and `allocation_completeness` at the end derives the required rows
# from the export list — a fifth in-place operation fails the suite until it is measured here or exempted
# with a reason. The bound is strict (`< bound`), so an exact zero is asserted as `< 1`.
@testset "allocation audit" begin
  Random.seed!(1)
  batch = randn(10_000)

  h = Histogram1D(-8.0, 8.0, 1000; order = 6)
  # batch observe! is the KA scatter kernel: a bounded launch cost, no O(n)/O(nbins) temporary.
  @alloc_below 4096 observe!(h, batch)
  # the scalar fold is a direct host update: exactly zero
  @alloc_below 1 observe!(h, 0.4)

  summary_a = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(summary_a, batch)
  summary_b = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(summary_b, batch)
  @alloc_below 1 merge!(summary_a, summary_b)

  # the cone operations (scale! and the weighted merge!) are in-place broadcasts — also allocation-free
  c = Histogram1D(-8.0, 8.0, 1000; order = 6, count_type = Float64); observe!(c, batch)
  f = Histogram1D(-8.0, 8.0, 1000; order = 6, count_type = Float64); observe!(f, batch)
  @alloc_below 1 scale!(c, 0.9)
  @alloc_below 1 merge!(c, f; weight = 0.5)

  # the Histogram2D hot path: the same bounded-launch contract for the pair scatter, zero for merge!
  xb = randn(10_000); yb = randn(10_000)
  j = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80)
  @alloc_below 4096 observe_pairs!(j, xb, yb)
  ja = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(ja, xb, yb)
  jb = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(jb, xb, yb)
  @alloc_below 1 merge!(ja, jb)
end


# The readouts are off the hot path, but a consumer calls them every step on a query grid
# (GBGCloudPhysics' closure), so their allocation must not scale as bins × queries. Each may allocate its
# output, the bin locations, and an O(degree) fit per query — never a per-bin temporary per query. A
# kernel passed through a helper without its type parameter is called by dynamic dispatch, ~70 bytes per
# bin per query (13 MiB here); a CPU summary copied by `Array` costs its whole storage per call. The
# budget is 3 KiB per query plus 16 bytes per bin.
@testset "readout allocation: no per-bin cost per query" begin
  Random.seed!(1)
  summary = Histogram1D(-8.0, 8.0, 1000); observe!(summary, randn(10_000))
  xb = randn(10_000)
  joint = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(joint, xb, xb .+ randn(10_000))
  xs = collect(range(-2.0, 2.0; length = 200))
  budget_1d = 3072 * length(xs) + 16 * length(summary.counts)
  budget_2d = 3072 * length(xs) + 16 * length(joint.counts)
  @alloc_below budget_1d density(summary, xs; h = 0.3)
  @alloc_below budget_1d density(summary, xs; h = 0.3, degree = 1)
  @alloc_below budget_1d density_standard_error(summary, xs; h = 0.3)
  @alloc_below budget_2d conditional_mean(joint, xs; h = 0.4)
  @alloc_below budget_2d conditional_density(joint, 0.3, xs; h_x = 0.4, h_y = 0.3)
  @alloc_below budget_2d conditional_expectation_standard_error(joint, cos, xs; h = 0.4)
end

# Completeness, derived from the export list rather than remembered: every exported in-place operation
# must have been exercised above, or be named here with a reason. None is exempt — the four `!`
# operations are the whole hot path, and each has a row.
GBGCodingTools.TestHarness.allocation_completeness(GBGStreamingStatistics, PACKAGE_MODULES)
