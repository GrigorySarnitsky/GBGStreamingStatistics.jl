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
    @test bandwidth(h; rule = :scott) ≈ bw / (4 / 3)^(1 / 5)

    truepdf(x) = exp(-x^2 / 2) / sqrt(2π)
    for x in -3.0:0.5:3.0
      @test density(h, x; h = bw) ≈ truepdf(x) atol = 3e-3
    end

    # normalization: the estimate integrates to 1 (trapezoid on a fine grid covering the support)
    xs = collect(-8.0:0.01:8.0)
    @test sum(density(h, xs; h = bw)) * 0.01 ≈ 1 atol = 1e-2
  end

  @testset "Ornstein–Uhlenbeck oracle (conditional / transition, §4.3–4.4)" begin
    # Forthcoming with the conditional estimators (Histogram2D + observe_pairs!): drift, diffusion,
    # and the transition density p_τ(y∣x) are known in closed form for the OU process.
    @test_skip conditional_mean_matches_ou_oracle()
    @test_skip transition_density_matches_ou_oracle()
  end

  @testset "consistency / rate" begin
    # Forthcoming: the integrated error of f̂ falls at the theoretical n^{-2/5} as n grows.
    @test_skip estimation_error_falls_as_n_to_the_minus_two_fifths()
  end
end
