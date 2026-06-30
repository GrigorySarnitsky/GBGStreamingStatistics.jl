using Test
using GBGStreamingStatistics
using Random
include("testoracles.jl"); using .TestOracles
include("testcontracts.jl"); using .TestContracts

@testset "ConditionalDensity (f_{Y|X}, transition density §4.3)" begin
  # Bivariate-Gaussian fixture: Y = ρX + √(1−ρ²)·ε ⇒ Y | X=x ∼ N(ρx, 1−ρ²) — an exact analytic
  # conditional density to check normalization and pointwise recovery against.
  Random.seed!(11)
  N = 2_000_000
  ρ = 0.6
  X = randn(N)
  Y = ρ .* X .+ sqrt(1 - ρ^2) .* randn(N)
  g = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  observe_pairs!(g, X, Y)
  hx = 0.3; hy = 0.2

  cond_pdf(y, x) = gaussian_pdf(y; mu = ρ * x, sigma = sqrt(1 - ρ^2))   # N(ρx, 1−ρ²)

  @testset "normalization: ∫ f̂(y|x) dy ≈ 1" begin
    ys = collect(-6.0:0.02:6.0)
    for x in (-1.0, 0.0, 1.0)
      fy = conditional_density(g, x, ys; hx = hx, hy = hy)
      @test sum(fy) * 0.02 ≈ 1 atol = 1e-2
    end
  end

  @testset "recovers N(ρx, 1−ρ²)" begin
    for x in (-0.5, 0.0, 0.5), y in (-0.5, 0.0, 0.5)
      @test conditional_density(g, x, y; hx = hx, hy = hy) ≈ cond_pdf(y, x) atol = 0.03
    end
  end

  @testset "the kernel option (parity): Epanechnikov normalizes and recovers" begin
    # The double-kernel now takes any symmetric ∫K=1 kernel; Epanechnikov still integrates to 1 (the
    # ∫K=1 normalization is what guarantees it) and recovers the conditional to a few percent.
    ys = collect(-6.0:0.02:6.0)
    for x in (-0.5, 0.0, 0.5)
      fy = conditional_density(g, x, ys; hx = hx, hy = hy, kernel = epanechnikov_kernel)
      @test sum(fy) * 0.02 ≈ 1 atol = 1e-2
    end
    for x in (-0.5, 0.0, 0.5), y in (-0.5, 0.0, 0.5)
      @test conditional_density(g, x, y; hx = hx, hy = hy, kernel = epanechnikov_kernel) ≈
            cond_pdf(y, x) atol = 0.05
      @test transition_density(g, x, y; hx = hx, hy = hy, kernel = epanechnikov_kernel) ≈
            conditional_density(g, x, y; hx = hx, hy = hy, kernel = epanechnikov_kernel)
    end
  end

  @testset "estimand-specific API: transition_density alias, NaN off the conditioning support, bandwidth errors" begin
    ys = [-0.5, 0.0, 0.5]
    fv = conditional_density(g, 0.0, ys; hx = hx, hy = hy)
    @test transition_density(g, 0.0, 0.0; hx = hx, hy = hy) ≈ fv[2]      # alias = conditional
    @test transition_density(g, 0.0, ys; hx = hx, hy = hy) == fv
    @test isnan(conditional_density(g, 100.0, 0.0; hx = hx, hy = hy))    # off-support conditioning x ⇒ NaN
    @test_throws ArgumentError conditional_density(g, 0.0, 0.0; hx = -1.0, hy = hy)
    @test_throws ArgumentError conditional_density(g, 0.0, 0.0; hx = hx, hy = 0.0)
  end

  # the universal query-point contract (scalar==vector over y, a functional of the merged summary),
  # applied by construction across a range of conditioning x — see test/testcontracts.jl
  half = N ÷ 2
  ga = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  gb = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  observe_pairs!(ga, view(X, 1:half), view(Y, 1:half))
  observe_pairs!(gb, view(X, half+1:N), view(Y, half+1:N))
  merge!(ga, gb)
  for x0 in (-0.5, 0.0, 0.5)
    check_point_estimator((s, y) -> conditional_density(s, x0, y; hx = hx, hy = hy), ga, g; points = [-0.3, 0.0, 0.3], rtol = 1e-9)
  end
end
