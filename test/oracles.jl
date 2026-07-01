using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random, Statistics
const H = Histograms
include("testoracles.jl"); using .TestOracles

@testset "analytic oracles (Gaussian, Ornstein–Uhlenbeck)" begin
  @testset "density recovers the Gaussian pdf (§4.1)" begin
    Random.seed!(2024)
    data = randn(2_000_000)
    h = Histogram1D(-8.0, 8.0, 2000; order = 4)
    observe!(h, data)

    bw = bandwidth(h)
    @test bw > 0
    @test bw ≈ (4 / 3)^(1 / 5) * GBGStreamingStatistics.Expectations._sample_std(h) * float(H._total(h))^(-1 / 5)
    @test bandwidth(h; kernel = epanechnikov_kernel) ≈ 2.2138 * bw  rtol = 2e-3   # kernel-aware: δ(Epan)/δ(Gauss)

    for x in -3.0:0.5:3.0
      @test density(h, x; h = bw) ≈ gaussian_pdf(x) atol = 3e-3
    end

    # normalization: the estimate integrates to 1 (trapezoid on a fine grid covering the support)
    xs = collect(-8.0:0.01:8.0)
    @test sum(density(h, xs; h = bw)) * 0.01 ≈ 1 atol = 1e-2
  end

  @testset "Ornstein–Uhlenbeck oracle: transition density (§4.3)" begin
    # OU process  dX = −θ X dt + σ dW : stationary N(0, s²) with s² = σ²/(2θ); the lag-τ transition is
    # exactly  X_{t+τ} | X_t = x ∼ N(x e^{−θτ}, s²(1 − e^{−2θτ})), so exact stationary lagged pairs are
    # sampled directly (no Euler discretization). (The §4.4 Kramers–Moyal drift/diffusion reading of the
    # same pairs is archived — see archive/.)
    Random.seed!(2025)
    θ = 0.8; σ = 1.2; τ = 0.5
    s = ou_stationary_std(theta = θ, sigma = σ)
    n = 3_000_000
    x0 = s .* randn(n)
    xτ = x0 .* exp(-θ * τ) .+ (s * sqrt(1 - exp(-2θ * τ))) .* randn(n)
    g = Histogram2D(-5s, 5s, 300, -5s, 5s, 300; lag = 1)
    observe_pairs!(g, x0, xτ)

    # transition density (§4.3): p_τ(y∣x) = N(x e^{−θτ}, s²(1−e^{−2θτ})), from the OU pairs. The
    # double-kernel estimate carries the usual O(h²) KDE bias, worst at the mode; atol reflects it.
    ptrans(y, x) = ou_transition(y, x; theta = θ, sigma = σ, tau = τ)
    for x in (-0.5, 0.0, 0.5)
      μx = x * exp(-θ * τ)
      for y in (μx - 0.4, μx, μx + 0.4)
        @test transition_density(g, x, y; hx = 0.25, hy = 0.15) ≈ ptrans(y, x) atol = 0.04
      end
    end
    ys = collect(-5s:0.02:5s)
    @test sum(transition_density(g, 0.0, ys; hx = 0.25, hy = 0.15)) * 0.02 ≈ 1 atol = 1e-2
  end
  # the density convergence-rate test (ISE ∝ n^{-2p/(2p+1)} per degree) and the representative
  # analytic-PDF accuracy suite live in `test/density_convergence.jl`.
end
