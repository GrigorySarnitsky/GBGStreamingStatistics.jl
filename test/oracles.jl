using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random, Statistics
const H = Histograms

@testset "analytic oracles (Gaussian, Ornstein–Uhlenbeck)" begin
  @testset "density recovers the Gaussian pdf (§4.1)" begin
    Random.seed!(2024)
    data = randn(2_000_000)
    h = Histogram1D(-8.0, 8.0, 2000; order = 4)
    observe!(h, data)

    bw = bandwidth(h)
    @test bw > 0
    @test bw ≈ (4 / 3)^(1 / 5) * H._sample_std(h) * float(H._total(h))^(-1 / 5)
    @test bandwidth(h; kernel = epanechnikov_kernel) ≈ 2.2138 * bw  rtol = 2e-3   # kernel-aware: δ(Epan)/δ(Gauss)

    truepdf(x) = exp(-x^2 / 2) / sqrt(2π)
    for x in -3.0:0.5:3.0
      @test density(h, x; h = bw) ≈ truepdf(x) atol = 3e-3
    end

    # normalization: the estimate integrates to 1 (trapezoid on a fine grid covering the support)
    xs = collect(-8.0:0.01:8.0)
    @test sum(density(h, xs; h = bw)) * 0.01 ≈ 1 atol = 1e-2
  end

  @testset "Ornstein–Uhlenbeck oracle: drift, diffusion & transition density (§4.3–4.4)" begin
    # OU process  dX = −θ X dt + σ dW : stationary N(0, s²) with s² = σ²/(2θ); the lag-τ transition is
    # exactly  X_{t+τ} | X_t = x ∼ N(x e^{−θτ}, s²(1 − e^{−2θτ})), so exact stationary lagged pairs are
    # sampled directly (no Euler discretization). The finite-τ Kramers–Moyal coefficients are known:
    #   drift     D⁽¹⁾(x,τ) = x (e^{−θτ} − 1)/τ                          → −θx    as τ → 0,
    #   diffusion D⁽²⁾(x,τ) = [s²(1 − e^{−2θτ}) + x²(e^{−θτ} − 1)²]/(2τ)  → σ²/2  as τ → 0.
    Random.seed!(2025)
    θ = 0.8; σ = 1.2; τ = 0.5
    s = sqrt(σ^2 / (2θ))
    n = 3_000_000
    x0 = s .* randn(n)
    xτ = x0 .* exp(-θ * τ) .+ (s * sqrt(1 - exp(-2θ * τ))) .* randn(n)
    g = Histogram2D(-5s, 5s, 300, -5s, 5s, 300; lag = 1)
    observe_pairs!(g, x0, xτ)

    D1(x) = x * (exp(-θ * τ) - 1) / τ
    D2(x) = (s^2 * (1 - exp(-2θ * τ)) + x^2 * (exp(-θ * τ) - 1)^2) / (2τ)
    hh = 0.2
    for x in (-1.0, -0.5, 0.0, 0.5, 1.0)
      # drift is essentially unbiased — the OU conditional mean E[Y|x] = x e^{−θτ} is linear, so the
      # local-linear fit has no smoothing bias.
      @test drift(g, x; h = hh, τ = τ) ≈ D1(x) atol = 0.02
      # diffusion carries the local-linear smoothing bias of the curved E[Y²|x] (quadratic in x), as
      # for conditional_variance — a few percent at this bandwidth (Fan–Yao refinement forthcoming).
      @test diffusion(g, x; h = hh, τ = τ) ≈ D2(x) rtol = 0.08
    end
    # the increment moments line up: M_1 = τ D⁽¹⁾, M_2 = 2τ D⁽²⁾
    @test increment_moment(g, 0.5, 1; h = hh) ≈ τ * D1(0.5)  atol = 0.02
    @test increment_moment(g, 0.5, 2; h = hh) ≈ 2τ * D2(0.5) rtol = 0.08

    # transition density (§4.3): p_τ(y∣x) = N(x e^{−θτ}, s²(1−e^{−2θτ})), from the same OU pairs. The
    # double-kernel estimate carries the usual O(h²) KDE bias, worst at the mode; atol reflects it.
    vt = s^2 * (1 - exp(-2θ * τ))
    ptrans(y, x) = exp(-(y - x * exp(-θ * τ))^2 / (2vt)) / sqrt(2π * vt)
    for x in (-0.5, 0.0, 0.5)
      μx = x * exp(-θ * τ)
      for y in (μx - 0.4, μx, μx + 0.4)
        @test transition_density(g, x, y; hx = 0.25, hy = 0.15) ≈ ptrans(y, x) atol = 0.04
      end
    end
    ys = collect(-5s:0.02:5s)
    @test sum(transition_density(g, 0.0, ys; hx = 0.25, hy = 0.15)) * 0.02 ≈ 1 atol = 1e-2
  end

  @testset "consistency / rate" begin
    # Forthcoming: the integrated error of f̂ falls at the theoretical n^{-2/5} as n grows.
    @test_skip estimation_error_falls_as_n_to_the_minus_two_fifths()
  end
end
