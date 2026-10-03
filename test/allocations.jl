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
#
# THE TRUE ZERO IS PINNED BY A PLAIN `@allocated`, as the rule asks where a bound covers something else.
# `@alloc_below` wraps its measurement in `@test (@allocated …) < limit`, and on Julia 1.12 that wrapping
# itself reads 48–208 bytes for the `merge!` calls (measured 2026-10-03: a plain `@allocated merge!(a, b)`
# is 0 in global, function, `try` and `@testset` scope alike; only the macro's context is not). So each
# `merge!` row carries a small bound for the audit, and the exact-zero assertion stands beside it; the
# macro's overhead is reported to GBGCodingTools.
@testset "allocation audit" begin
  Random.seed!(1)
  batch = randn(10_000)

  h = Histogram1D(-8.0, 8.0, 1000; order = 6)
  # batch observe! is the KA scatter kernel: a bounded launch cost, no O(n)/O(nbins) temporary.
  @alloc_below 4096 observe!(h, batch)
  # the scalar fold is a direct host update: exactly zero
  @alloc_below 1 observe!(h, 0.4)

  a = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(a, batch)
  b = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(b, batch)
  @test (@allocated merge!(a, b)) == 0                  # the true zero (merge! is warm from the line above)
  @alloc_below 256 merge!(a, b)                         # the audit row; the bound is the harness's own overhead

  # the cone operations (scale! and the weighted merge!) are in-place broadcasts — also allocation-free
  c = Histogram1D(-8.0, 8.0, 1000; order = 6, count = Float64); observe!(c, batch)
  f = Histogram1D(-8.0, 8.0, 1000; order = 6, count = Float64); observe!(f, batch)
  @alloc_below 1 scale!(c, 0.9)
  merge!(c, f; weight = 0.5)                            # warm up the weighted method
  @test (@allocated merge!(c, f; weight = 0.5)) == 0    # the true zero
  @alloc_below 256 merge!(c, f; weight = 0.5)

  # the Histogram2D hot path: the same bounded-launch contract for the pair scatter, zero for merge!
  xb = randn(10_000); yb = randn(10_000)
  j = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80)
  @alloc_below 4096 observe_pairs!(j, xb, yb)
  ja = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(ja, xb, yb)
  jb = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(jb, xb, yb)
  merge!(ja, jb)                                        # warm up
  @test (@allocated merge!(ja, jb)) == 0                # the true zero
  @alloc_below 256 merge!(ja, jb)
end

# Completeness, derived from the export list rather than remembered: every exported in-place operation
# must have been exercised above, or be named here with a reason. None is exempt — the four `!`
# operations are the whole hot path, and each has a row.
GBGCodingTools.TestHarness.allocation_completeness(GBGStreamingStatistics, PACKAGE_MODULES)
