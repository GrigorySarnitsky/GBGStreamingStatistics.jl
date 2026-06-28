using Test
using GBGStreamingStatistics
using Random

# The hot path is allocation-free (invariant 4): on a preallocated histogram, observe!/merge! allocate
# no bytes (measured on the CPU via @allocated, after a warm-up call to force compilation). This is the
# contract that keeps a long GPU run from wedging once the device kernel lands.
@testset "allocation audit (@allocated)" begin
  Random.seed!(1)
  batch = randn(10_000)

  h = Histogram1D(-8.0, 8.0, 1000; order = 6)
  observe!(h, batch)                                   # warm up (compile)
  # batch observe! is the KA scatter kernel: a bounded launch cost, no O(n)/O(nbins) temporary.
  @test (@allocated observe!(h, batch)) < 4096

  observe!(h, 0.3)                                      # warm up the scalar (direct host) method
  @test (@allocated observe!(h, 0.4)) == 0

  a = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(a, batch)
  b = Histogram1D(-8.0, 8.0, 1000; order = 6); observe!(b, batch)
  merge!(a, b)                                          # warm up
  @test (@allocated merge!(a, b)) == 0

  # the Histogram2D hot path (observe_pairs! / merge!) has the same bounded-allocation contract; it is
  # asserted alongside the 2D summary in test/histograms2d.jl.
end
