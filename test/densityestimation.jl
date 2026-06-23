using Test
using GBGStreamingStatistics

@testset "DensityEstimation (density)" begin
  K(u) = exp(-u^2 / 2) / sqrt(2π)

  # One point binned at a known bin center → the binned KDE is a single Gaussian bump at that center
  # (the WARPing approximation evaluates at bin centers, not the raw sample).
  h = Histogram1D(-1.0, 1.0, 2)        # bins [-1,0) center −0.5, [0,1) center 0.5
  observe!(h, 0.4)                      # lands in bin 2 (center 0.5)
  bw = 0.2
  @test density(h, 0.5; h = bw) ≈ K(0.0) / bw
  @test density(h, 0.3; h = bw) ≈ K((0.3 - 0.5) / bw) / bw

  # an empty histogram has zero density everywhere
  @test density(Histogram1D(-1.0, 1.0, 4), 0.0; h = 0.1) == 0.0

  # a vector of query points returns a vector
  ys = density(h, [0.3, 0.5, 0.7]; h = bw)
  @test ys isa AbstractVector
  @test length(ys) == 3
  @test ys[2] ≈ density(h, 0.5; h = bw)

  # the bandwidth must be positive
  @test_throws ArgumentError density(h, 0.0; h = 0.0)
  @test_throws ArgumentError density(h, 0.0; h = -0.5)
end
