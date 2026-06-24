using Test
using GBGStreamingStatistics
using Random

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

  cond_pdf(y, x) = exp(-(y - ρ * x)^2 / (2 * (1 - ρ^2))) / sqrt(2π * (1 - ρ^2))   # N(ρx, 1−ρ²)

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

  @testset "API: scalar/vector forms, transition_density alias, NaN off-support" begin
    ys = [-0.5, 0.0, 0.5]
    fv = conditional_density(g, 0.0, ys; hx = hx, hy = hy)
    @test fv isa AbstractVector && length(fv) == 3
    @test fv[2] ≈ conditional_density(g, 0.0, 0.0; hx = hx, hy = hy)     # vector vs scalar
    @test transition_density(g, 0.0, 0.0; hx = hx, hy = hy) ≈ fv[2]      # alias = conditional
    @test transition_density(g, 0.0, ys; hx = hx, hy = hy) == fv
    @test isnan(conditional_density(g, 100.0, 0.0; hx = hx, hy = hy))    # no local data ⇒ NaN
    @test_throws ArgumentError conditional_density(g, 0.0, 0.0; hx = -1.0, hy = hy)
    @test_throws ArgumentError conditional_density(g, 0.0, 0.0; hx = hx, hy = 0.0)
  end

  @testset "mergeability: the density is a functional of the merged summary" begin
    half = N ÷ 2
    ga = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
    gb = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
    observe_pairs!(ga, view(X, 1:half), view(Y, 1:half))
    observe_pairs!(gb, view(X, half+1:N), view(Y, half+1:N))
    merge!(ga, gb)
    for (x, y) in ((-0.5, 0.0), (0.0, 0.0), (0.5, 0.3))
      @test conditional_density(ga, x, y; hx = hx, hy = hy) ≈
            conditional_density(g, x, y; hx = hx, hy = hy) rtol = 1e-10
    end
  end
end
