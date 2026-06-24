using Test
using GBGStreamingStatistics
using Random

@testset "MarkovInference (increment moments, drift, diffusion §4.4)" begin
  # A bivariate-Gaussian fixture (Y = bX + noise standing in for the lagged response) suffices to
  # check the API contract, the algebraic identities, and mergeability; the analytic correctness of
  # drift/diffusion against the Ornstein–Uhlenbeck closed forms lives in test/oracles.jl.
  Random.seed!(7)
  N = 1_000_000
  τ = 2
  X = randn(N)
  Y = 0.5 .* X .+ sqrt(0.75) .* randn(N)
  g = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; lag = τ)
  observe_pairs!(g, X, Y)
  hh = 0.4

  @testset "M_1, M_2 and the variance decomposition" begin
    x = 0.5
    M1 = increment_moment(g, x, 1; h = hh)
    M2 = increment_moment(g, x, 2; h = hh)
    # M_1(x) = E[Y|x] − x
    @test M1 ≈ conditional_mean(g, x; h = hh) - x atol = 1e-9
    # M_2(x) = E[(Y−x)²|x] = Var(Y|x) + M_1²   (the variance decomposition)
    @test M2 ≈ conditional_variance(g, x; h = hh) + M1^2 rtol = 1e-6
    ms = increment_moment(g, [-0.5, 0.0, 0.5], 1; h = hh)
    @test ms isa AbstractVector && length(ms) == 3
  end

  @testset "Kramers–Moyal scaling: drift = D⁽¹⁾, diffusion = D⁽²⁾" begin
    x = 0.3
    @test drift(g, x; h = hh)     ≈ increment_moment(g, x, 1; h = hh) / τ     rtol = 1e-12
    @test diffusion(g, x; h = hh) ≈ increment_moment(g, x, 2; h = hh) / (2τ)  rtol = 1e-12
    @test drift(g, x; h = hh)     ≈ kramers_moyal(g, x, 1; h = hh)            rtol = 1e-12
    @test diffusion(g, x; h = hh) ≈ kramers_moyal(g, x, 2; h = hh)            rtol = 1e-12
    # the lag time τ enters as 1/(n!·τ): doubling τ halves D⁽¹⁾
    @test drift(g, x; h = hh, τ = 2τ) ≈ drift(g, x; h = hh) / 2 rtol = 1e-12
    ds = drift(g, [-0.5, 0.0, 0.5]; h = hh)
    @test ds isa AbstractVector && length(ds) == 3
  end

  @testset "coverage gaps and argument checks" begin
    @test_throws ArgumentError increment_moment(g, 0.0, 3; h = hh)    # n ≥ 3 needs order > 2
    @test_throws ArgumentError increment_moment(g, 0.0, 0; h = hh)    # n < 1
    @test_throws ArgumentError increment_moment(g, 0.0, 1; h = -1.0)  # bad bandwidth
    g0 = Histogram2D(-6.0, 6.0, 60, -6.0, 6.0, 60)                    # lag = 0 (static joint)
    observe_pairs!(g0, X, Y)
    @test_throws ArgumentError drift(g0, 0.0; h = hh)                 # no positive lag time
    @test drift(g0, 0.0; h = hh, τ = 1.0) isa Real                    # explicit τ is fine
  end

  @testset "mergeability: drift/diffusion are functionals of the merged summary" begin
    half = N ÷ 2
    ga = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; lag = τ)
    gb = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; lag = τ)
    observe_pairs!(ga, view(X, 1:half), view(Y, 1:half))
    observe_pairs!(gb, view(X, half+1:N), view(Y, half+1:N))
    merge!(ga, gb)
    for x in (-0.5, 0.0, 0.5)
      @test drift(ga, x; h = hh)     ≈ drift(g, x; h = hh)     rtol = 1e-8
      @test diffusion(ga, x; h = hh) ≈ diffusion(g, x; h = hh) rtol = 1e-8
    end
  end
end
