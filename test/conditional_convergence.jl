using Test
using GBGStreamingStatistics
using Random, Statistics
include("testoracles.jl"); using .TestOracles

# The nonparametric convergence rates of the conditional estimators, mirroring the density ISE-slope test
# (test/density_convergence.jl): accumulate an analytic oracle at growing n, hold the bandwidth at the
# rate-optimal scaling, fit the log–log slope of the integrated squared error, and check it against the
# theoretical exponent. Two genuinely distinct rates:
#   conditional_mean   — the local-linear L²-projection m(x) = E[Y|X=x]: pointwise RMSE ∝ n^{-2/5}, so
#                        ISE ∝ n^{-4/5} (slope −0.80), the same nonparametric-regression rate as a 1-D KDE.
#   conditional_density — the double-kernel f(y|x) is a 2-D smoothing problem (d = 2): RMISE ∝ n^{-1/3},
#                        so ISE ∝ n^{-2/3} (slope −0.667), the slower 2-D rate.
# drift / diffusion inherit conditional_mean's rate exactly — drift = (m̂(x) − x)/τ is an affine functional
# of the same local-linear fit, and diffusion a functional of the first two conditional moments at the same
# rate — so they need no separate slope test (docs/TESTING.md coverage matrix). The tolerances and seed
# blocks are audited across seeds (the threshold sits above the worst block with margin), as in
# docs/TESTING.md "Tolerances are seed-robust".

loglog_slope(xs, ys) = let lx = log.(xs), ly = log.(ys), mx = mean(log.(xs)), my = mean(log.(ys))
  sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
end


@testset "conditional_mean: ISE ∝ n^(−4/5) (the n^(−2/5) L²-projection rate, §4.2)" begin
  # Curved-mean oracle: X ~ N(0,1), m(x) = sin(x) (m'' ≢ 0, so the local-linear fit carries a genuine
  # O(h²) bias — a linear mean would be fit exactly and would NOT exhibit the nonparametric rate, the
  # curvature is the point), Y = sin(X) + σε. At the rate-optimal bandwidth h ∝ n^{-1/5} the pointwise
  # MSE balances bias² ~ h⁴ and variance ~ 1/(nh), both n^{-4/5}; the ISE over the interior inherits it.
  σ = 0.35
  c = 1.25                                       # h = c·n^{-1/5}, the degree-1 rate-optimal scaling
  ns = (30_000, 80_000, 200_000, 500_000, 1_200_000)
  grid = collect(-2.0:0.04:2.0); truth = sin.(grid)
  function mean_ise(n; seed)
    Random.seed!(seed)
    X = randn(n); Y = sin.(X) .+ σ .* randn(n)
    g = Histogram2D(-6.0, 6.0, 600, -4.0, 4.0, 160); observe_pairs!(g, X, Y)
    return ise(conditional_mean(g, grid; h = c * n^(-1 / 5)), truth, 0.04)
  end
  ises = [mean(mean_ise(n; seed = 100 + s) for s in 1:4) for n in ns]
  slope = loglog_slope(collect(ns), ises)
  @test slope ≈ -0.80 atol = 0.12      # ISE ∝ n^{-4/5}; RMSE ∝ n^{-2/5}, the L²-projection rate
  @test ises[end] < 0.2 * ises[1]      # and a decisive drop across the n-decades
end


@testset "conditional_density: ISE ∝ n^(−2/3) (the 2-D double-kernel rate, §4.3)" begin
  # Bivariate-Gaussian oracle: Y | X=x ∼ N(ρx, 1−ρ²), smooth and curved in both x and y, so the
  # double-kernel estimate carries the 2-D O(hx²+hy²) bias. At the 2-D rate-optimal bandwidth
  # hx = hy ∝ n^{-1/(d+4)} = n^{-1/6} the MISE balances bias² ~ h⁴ and variance ~ 1/(n·hx·hy), both
  # n^{-2/3}. The ISE is measured in y at the data-densest conditioning point x₀ = 0 (f_X maximal ⇒
  # cleanest slope). The n-range starts higher than the 1-D density test: the slower 2-D rate reaches its
  # asymptotic regime later, and below ~10^5 the measured slope runs shallow (≈ −0.6), the conditional
  # analog of the finite-n muddying the density rate test records for non-Gaussians.
  ρ = 0.6
  sd = sqrt(1 - ρ^2)
  c = 1.5                                        # hx = hy = c·n^{-1/6}, the 2-D rate-optimal scaling
  ns = (130_000, 330_000, 800_000, 1_500_000, 2_800_000)
  ys = collect(-3.0:0.03:3.0); truth = gaussian_pdf.(ys; mu = 0.0, sigma = sd)
  function cd_ise(n; seed)
    Random.seed!(seed)
    X = randn(n); Y = ρ .* X .+ sd .* randn(n)
    g = Histogram2D(-6.0, 6.0, 400, -6.0, 6.0, 500; lag = 1); observe_pairs!(g, X, Y)
    hb = c * n^(-1 / 6)
    return ise(conditional_density(g, 0.0, ys; hx = hb, hy = hb), truth, 0.03)
  end
  ises = [mean(cd_ise(n; seed = 200 + s) for s in 1:4) for n in ns]
  slope = loglog_slope(collect(ns), ises)
  @test slope ≈ -0.667 atol = 0.10     # MISE ∝ n^{-2/3}; RMISE ∝ n^{-1/3}, the 2-D double-kernel rate
  @test ises[end] < 0.4 * ises[1]      # and a clear drop (the slower 2-D rate over a ~20× n-range)
end
