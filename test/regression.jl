using Test
using GBGStreamingStatistics
using Random
include("testcontracts.jl"); using .TestContracts

@testset "Regression (conditional_mean / conditional_variance)" begin
  # Bivariate-Gaussian oracle: Y = ρX + √(1−ρ²)·ε ⇒ E[Y|X=x] = ρx (linear), Var(Y|X=x) = 1−ρ² (const).
  Random.seed!(2024)
  N = 2_000_000
  ρ = 0.6
  X = randn(N)
  Y = ρ .* X .+ sqrt(1 - ρ^2) .* randn(N)
  hist = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120)
  observe_pairs!(hist, X, Y)

  @testset "conditional_mean recovers ρx (local-linear)" begin   # estimand-specific: the analytic oracle
    for x in -2.0:0.5:2.0
      @test conditional_mean(hist, x; h = 0.4) ≈ ρ * x atol = 0.02
    end
    @test_throws ArgumentError conditional_mean(hist, 0.0; h = -1.0)
  end

  @testset "conditional_variance recovers 1−ρ² (const)" begin    # estimand-specific: the analytic oracle
    # The simple Ê[Y²|x]−m̂² estimator carries the local-linear smoothing bias of the curved
    # E[Y²|x]; at a modest bandwidth it is within a few percent (Fan–Yao refinement forthcoming).
    for x in (-1.0, 0.0, 1.0)
      @test conditional_variance(hist, x; h = 0.2) ≈ 1 - ρ^2 rtol = 0.1
    end
    @test conditional_variance(hist, 0.0; h = 0.2) > 0
  end

  @testset "degree and kernel knobs (parity with density)" begin   # the new local_fit-backed knobs
    # The bivariate-Gaussian mean is linear, so local-quadratic (degree 2) recovers it, and a second
    # kernel (Epanechnikov) agrees with the Gaussian default. Nadaraya–Watson (degree 0) carries the
    # design-density bias 2 m' f'/f away from the centre, so it is checked only where that vanishes
    # (x = 0); degree outside {0, 1, 2} is rejected.
    for x in (-1.0, 0.0, 1.0)
      @test conditional_mean(hist, x; h = 0.5, degree = 2) ≈ ρ * x atol = 0.03
      @test conditional_mean(hist, x; h = 0.4, kernel = epanechnikov_kernel) ≈ ρ * x atol = 0.03
      @test isfinite(conditional_mean(hist, x; h = 0.4, degree = 0))
    end
    @test conditional_mean(hist, 0.0; h = 0.4, degree = 0) ≈ 0.0 atol = 0.02   # NW unbiased at the centre
    @test_throws ArgumentError conditional_mean(hist, 0.0; h = 0.4, degree = 3)
    @test_throws ArgumentError conditional_variance(hist, 0.0; h = 0.2, degree = 3)
  end

  # the universal query-point contract (scalar==vector, a functional of the merged summary, NaN off the
  # data's support), applied by construction to both conditional moments — see test/testcontracts.jl
  half = N ÷ 2
  a = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120)
  b = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120)
  observe_pairs!(a, view(X, 1:half), view(Y, 1:half))
  observe_pairs!(b, view(X, half+1:N), view(Y, half+1:N))
  merge!(a, b)
  check_point_estimator((s, x) -> conditional_mean(s, x; h = 0.4),     a, hist; points = [-1.0, 0.0, 1.0], outside = 100.0)
  check_point_estimator((s, x) -> conditional_variance(s, x; h = 0.2), a, hist; points = [-1.0, 0.0, 1.0], outside = 100.0)
end
