using Test
using GBGStreamingStatistics
using Random

@testset "Regression (conditional_mean / conditional_variance)" begin
  # Bivariate-Gaussian oracle: Y = ρX + √(1−ρ²)·ε ⇒ E[Y|X=x] = ρx (linear), Var(Y|X=x) = 1−ρ² (const).
  Random.seed!(2024)
  N = 2_000_000
  ρ = 0.6
  X = randn(N)
  Y = ρ .* X .+ sqrt(1 - ρ^2) .* randn(N)
  hist = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120)
  observe_pairs!(hist, X, Y)

  @testset "conditional_mean recovers ρx (local-linear)" begin
    for x in -2.0:0.5:2.0
      @test conditional_mean(hist, x; h = 0.4) ≈ ρ * x atol = 0.02
    end
    ms = conditional_mean(hist, collect(-2.0:1.0:2.0); h = 0.4)
    @test ms isa AbstractVector
    @test length(ms) == 5
    @test ms[3] ≈ 0 atol = 0.02                       # x = 0
    @test isnan(conditional_mean(hist, 100.0; h = 0.4))   # no local data ⇒ NaN
    @test_throws ArgumentError conditional_mean(hist, 0.0; h = -1.0)
  end

  @testset "conditional_variance recovers 1−ρ² (const)" begin
    # The simple Ê[Y²|x]−m̂² estimator carries the local-linear smoothing bias of the curved
    # E[Y²|x]; at a modest bandwidth it is within a few percent (Fan–Yao refinement forthcoming).
    for x in (-1.0, 0.0, 1.0)
      @test conditional_variance(hist, x; h = 0.2) ≈ 1 - ρ^2 rtol = 0.1
    end
    @test conditional_variance(hist, 0.0; h = 0.2) > 0
    vs = conditional_variance(hist, [-1.0, 0.0, 1.0]; h = 0.2)
    @test length(vs) == 3
  end
end
