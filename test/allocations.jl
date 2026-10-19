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

# Completeness, derived from the export list rather than remembered: every exported in-place operation
# must have been exercised above, or be named here with a reason. None is exempt — the four `!`
# operations are the whole hot path, and each has a row.
GBGCodingTools.TestHarness.allocation_completeness(GBGStreamingStatistics, PACKAGE_MODULES)
