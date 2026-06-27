using Test
using GBGStreamingStatistics
const H = GBGStreamingStatistics.Histograms

@testset "Bandwidth (bandwidth selection)" begin
  # A known sample: σ̂ recovered exactly from the binned moments.
  h = Histogram1D(0.0, 10.0, 100)
  observe!(h, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
  n = H._total(h)
  σ = H._sample_std(h)
  @test σ ≈ sqrt(60 / 9)                                # exact recovery from the binned moments

  # rule = :normal (default), Gaussian kernel, degree 0/1 → Silverman (4/3)^(1/5)·σ̂·n^(-1/5).
  @test bandwidth(h) ≈ (4 / 3)^(1 / 5) * σ * float(n)^(-1 / 5)  rtol = 1e-3
  @test bandwidth(h; degree = 1) ≈ bandwidth(h)  rtol = 1e-6           # deg 0,1 share the 2nd-order rule

  # the kernel enters only through the canonical factor δ(K): the Epanechnikov h is 2.214× the
  # Gaussian, i.e. ≈ 2.345·σ̂·n^(-1/5) (= statsmodels' normal_reference_constant).
  @test bandwidth(h; kernel = epanechnikov_kernel) ≈ 2.3449 * σ * float(n)^(-1 / 5)  rtol = 2e-3
  @test bandwidth(h; kernel = epanechnikov_kernel) / bandwidth(h) ≈ 2.2138  rtol = 2e-3

  # degree 2 is a higher-order (O(h⁴)) fit → the n^(-1/9) rate and a wider, more slowly shrinking h.
  # Matching the closed form tests both the constant and the -1/9 exponent.
  @test bandwidth(h; degree = 2) ≈ 1.0834 * σ * float(n)^(-1 / 9)  rtol = 3e-3
  @test bandwidth(h; degree = 2) > bandwidth(h; degree = 0)

  # errors
  @test_throws ArgumentError bandwidth(h; rule = :nope)
  @test_throws ArgumentError bandwidth(h; degree = -1)
  @test_throws ArgumentError bandwidth(h; degree = 5)                  # coverage gap: degree ≤ 3
  @test_throws ArgumentError bandwidth(Histogram1D(0.0, 1.0, 10))      # n = 0 < 2
  h1 = Histogram1D(0.0, 1.0, 10); observe!(h1, 0.5)
  @test_throws ArgumentError bandwidth(h1)                             # n = 1 < 2
end
