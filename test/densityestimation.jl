using Test
using GBGStreamingStatistics
using Random

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

@testset "kernels: Gaussian + Epanechnikov" begin
  # the two shipped kernels are public
  @test gaussian_kernel(0.0) ≈ 1 / sqrt(2π)
  @test epanechnikov_kernel(0.0) == 0.75                  # 3/4 at the mode
  @test epanechnikov_kernel(0.5) == 0.75 * (1 - 0.25)
  @test epanechnikov_kernel(1.0) == 0.0                   # compact support [-1, 1]
  @test epanechnikov_kernel(2.0) == 0.0

  h = Histogram1D(-1.0, 1.0, 2); observe!(h, 0.4)         # one sample, bin 2 (center 0.5)
  bw = 0.6
  @test density(h, 0.5; h = bw, kernel = epanechnikov_kernel) ≈ epanechnikov_kernel(0.0) / bw
  @test density(h, 0.5; h = bw, kernel = gaussian_kernel) ≈ density(h, 0.5; h = bw)   # explicit == default
  @test density(h, -0.9; h = 0.3, kernel = epanechnikov_kernel) == 0.0                # compact support → 0 far away

  # both kernels integrate the recovered density to ≈ 1 over the range
  Random.seed!(7)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(200_000))
  xs = range(-8.0, 8.0; length = 3201); dx = step(xs)
  for kern in (gaussian_kernel, epanechnikov_kernel)
    @test sum(density(hg, x; h = 0.3, kernel = kern) for x in xs) * dx ≈ 1 atol = 2e-2
  end
end
