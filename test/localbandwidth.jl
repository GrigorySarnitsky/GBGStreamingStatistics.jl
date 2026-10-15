using Test
using GBGStreamingStatistics
using Random

@testset "Local bandwidth h⋆(x) (response-aware plug-in, §4.5)" begin
  # Quadratic-mean Gaussian oracle: X ~ N(0,1), Y = X² + σ·ε. Then m(x) = x² (so m'' = 2, recovered
  # exactly by a degree-2 fit), σ²(x) = σ² (constant), f_X(x) = φ(x). The local AMISE bandwidth is then
  # known in closed form: h⋆(x) = [R(K)·σ² / (μ₂(K)²·4·φ(x)·n)]^{1/5}. Near the centre all three pilots
  # are clean; far from it the plug-in σ̂² carries the curved-E[Y²|x] bias (the Fan–Yao refinement).
  Random.seed!(11)
  n = 4_000_000
  σ = 0.5
  X = randn(n)
  Y = X .^ 2 .+ σ .* randn(n)
  hist = Histogram2D(-6.0, 6.0, 240, -8.0, 40.0, 240)
  observe_pairs!(hist, X, Y)

  RK = 1 / (2 * sqrt(π))                 # R(gaussian) = 1/(2√π)
  μ2 = 1.0                               # μ₂(gaussian) = 1
  φ(x) = exp(-x^2 / 2) / sqrt(2π)
  hstar_true(x) = (RK * σ^2 / (μ2^2 * 4.0 * φ(x) * n))^(1 / 5)
  hp = 0.4                               # explicit pilot (m'' = 2 exact for the quadratic; deterministic)

  @testset "recovers the closed-form h⋆(x) near the centre" begin
    for x in (-0.3, 0.0, 0.3)
      @test bandwidth(hist, x; rule = :local, pilot_bandwidth = hp) ≈ hstar_true(x) rtol = 0.12
    end
  end

  @testset "positive and finite on support; NaN off it; rejected rules" begin
    for x in -1.5:0.5:1.5
      h = bandwidth(hist, x; rule = :local, pilot_bandwidth = hp)
      @test isfinite(h) && h > 0
    end
    @test isnan(bandwidth(hist, 100.0; rule = :local, pilot_bandwidth = hp))     # off the data's support
    @test_throws ArgumentError bandwidth(hist, 0.0; rule = :sj)         # only :local on a joint summary
    @test_throws ArgumentError bandwidth(hist, 0.0; rule = :local, pilot_bandwidth = -1.0)
  end

  @testset "scalar == vector" begin
    xs = [-0.3, 0.0, 0.3, 0.8]
    hs = bandwidth(hist, xs; rule = :local, pilot_bandwidth = hp)
    @test hs ≈ [bandwidth(hist, x; rule = :local, pilot_bandwidth = hp) for x in xs]
  end

  @testset "default pilot (curvature-rate rule of thumb) is sane" begin
    # With no pilot_bandwidth the σ̂_X·n^{-1/9} pilot is used; h⋆ should still be finite, positive, and within a
    # small factor of the explicit-pilot value near the centre.
    for x in (-0.3, 0.0, 0.3)
      h = bandwidth(hist, x; rule = :local)
      @test isfinite(h) && h > 0
      @test 0.3 < h / bandwidth(hist, x; rule = :local, pilot_bandwidth = hp) < 3.0
    end
  end

  @testset "mergeable: h⋆ is a functional of the merged summary" begin
    a = Histogram2D(-6.0, 6.0, 240, -8.0, 40.0, 240)
    b = Histogram2D(-6.0, 6.0, 240, -8.0, 40.0, 240)
    half = n ÷ 2
    observe_pairs!(a, view(X, 1:half), view(Y, 1:half))
    observe_pairs!(b, view(X, half + 1:n), view(Y, half + 1:n))
    merge!(a, b)
    for x in (-0.3, 0.0, 0.5)
      @test bandwidth(a, x; rule = :local, pilot_bandwidth = hp) ≈
            bandwidth(hist, x; rule = :local, pilot_bandwidth = hp)
    end
  end
end
