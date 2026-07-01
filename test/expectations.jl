# Tests for the §4.B power-sum quadrature readouts: `expectation` (E[g(X)]), its sampling
# `expectation_standard_error`, and `pushforward_density` (the density of Y = g(X)). Oracles: the Gaussian
# moment-generating / trigonometric expectations (E[cos X] = e^{-σ²/2}, etc.), polynomial exactness of the
# quadrature, and the lognormal / χ²₁ pushforwards (the latter non-monotone). Plus the moment unification
# (E[x^p] equals the recovered moments) and mergeability.

using Test
using GBGStreamingStatistics
using Random, Statistics

const LOWER, UPPER, NBINS = -8.0, 8.0, 4096


@testset "expectation: Gaussian oracles + moment unification" begin
  Random.seed!(20260701)
  data = randn(300_000)                       # X ~ N(0, 1)
  h = Histogram1D(LOWER, UPPER, NBINS)
  observe!(h, data)

  # Trigonometric / MGF oracles: E[cos X] = e^{-1/2}, E[sin X] = 0, E[e^X] = e^{1/2}.
  @test expectation(h, cos) ≈ exp(-0.5) atol = 0.01
  @test abs(expectation(h, sin)) < 0.01
  @test expectation(h, exp) ≈ exp(0.5) rtol = 0.02

  # Moment unification — the g = x^p case *is* the recovered moment (interpolation path ⇒ closed form).
  @test expectation(h, identity) ≈ mean(h) rtol = 1e-10
  @test expectation(h, x -> x^2) ≈ var(h; corrected = false) + mean(h)^2 rtol = 1e-10
  @test expectation(h, x -> (x - mean(h))^3) ≈ central_moment(h, 3) rtol = 1e-8 atol = 1e-10

  # Empty summary ⇒ NaN.
  @test isnan(expectation(Histogram1D(LOWER, UPPER, NBINS), cos))
end


@testset "expectation: polynomial exactness (degree ≤ order)" begin
  # For a polynomial g of degree ≤ order the quadrature is exact (the power-sums capture the within-bin
  # moments), so it reproduces the raw empirical mean (1/n) Σ g(Xᵢ) to rounding — independent of bin width.
  data = [-2.3, -0.7, 0.1, 1.5, 3.2, -4.4, 2.9, 0.05]
  h = Histogram1D(LOWER, UPPER, NBINS)
  observe!(h, data)
  g(x) = 1.0 + 2.0x - 0.5x^2 + 0.1x^3         # degree 3 ≤ order 4
  @test expectation(h, g) ≈ mean(g.(data)) rtol = 1e-10

  # A smooth non-polynomial g: the binning error vs the raw empirical mean is tiny at fine resolution.
  Random.seed!(7)
  d = randn(20_000)
  hf = Histogram1D(LOWER, UPPER, NBINS)
  observe!(hf, d)
  @test expectation(hf, cos) ≈ mean(cos.(d)) atol = 1e-6
end


@testset "expectation_standard_error" begin
  Random.seed!(99)
  data = randn(120_000)
  h = Histogram1D(LOWER, UPPER, NBINS)
  observe!(h, data)
  n = float(length(data))

  # SE of the empirical mean: g = identity ⇒ std(uncorrected)/√n exactly.
  @test expectation_standard_error(h, identity) ≈ std(h; corrected = false) / sqrt(n) rtol = 1e-10

  # SE of E[cos X]: √(Var(cos X)/n), matching the raw empirical variance of cos(Xᵢ).
  @test expectation_standard_error(h, cos) ≈ sqrt(var(cos.(data); corrected = false) / n) atol = 1e-6
  @test expectation_standard_error(h, cos) > 0

  @test isnan(expectation_standard_error(Histogram1D(LOWER, UPPER, NBINS), cos))
end


@testset "pushforward_density: Y = exp(X) ~ Lognormal(0, 1)" begin
  Random.seed!(2024)
  data = randn(400_000)
  h = Histogram1D(LOWER, UPPER, NBINS)
  observe!(h, data)
  lognormal_pdf(y) = 1 / (y * sqrt(2π)) * exp(-(log(y))^2 / 2)

  for y in (0.5, 1.0, 2.0, 3.0)
    @test pushforward_density(h, exp, y; h = 0.12) ≈ lognormal_pdf(y) rtol = 0.06
  end

  # Exact normalization ∫ f̂ dy = 1 (symmetric kernel), by the trapezoidal rule over the support.
  ygrid = range(0.02, 20.0; length = 2000)
  curve = pushforward_density(h, exp, collect(ygrid); h = 0.12)
  integral = sum(curve) * step(ygrid) - 0.5 * (curve[1] + curve[end]) * step(ygrid)
  @test integral ≈ 1.0 atol = 0.03

  # Scalar and vector forms agree.
  @test pushforward_density(h, exp, [0.5, 1.0, 2.0]; h = 0.12) ≈
        [pushforward_density(h, exp, y; h = 0.12) for y in (0.5, 1.0, 2.0)]

  @test_throws ArgumentError pushforward_density(h, exp, 1.0; h = -0.1)
end


@testset "pushforward_density: Y = X² ~ χ²₁ (non-monotone g)" begin
  # A non-monotone transform: both preimages ±√y contribute automatically — no g⁻¹, no manual Jacobian.
  Random.seed!(303)
  data = randn(500_000)
  h = Histogram1D(LOWER, UPPER, NBINS)
  observe!(h, data)
  chisq1_pdf(y) = 1 / sqrt(2π * y) * exp(-y / 2)

  for y in (0.5, 1.0, 2.0, 3.0)                # away from the g'=0 caustic at 0
    @test pushforward_density(h, x -> x^2, y; h = 0.10) ≈ chisq1_pdf(y) rtol = 0.10
  end
end


@testset "expectation: mergeability (accumulate-then-merge = accumulate-together)" begin
  Random.seed!(11)
  d1 = randn(50_000)
  d2 = randn(70_000)
  a = Histogram1D(LOWER, UPPER, NBINS); observe!(a, d1)
  b = Histogram1D(LOWER, UPPER, NBINS); observe!(b, d2)
  full = Histogram1D(LOWER, UPPER, NBINS); observe!(full, vcat(d1, d2))
  merge!(a, b)
  @test expectation(a, cos) ≈ expectation(full, cos) rtol = 1e-12
  @test expectation(a, exp) ≈ expectation(full, exp) rtol = 1e-12
  @test pushforward_density(a, exp, 1.3; h = 0.15) ≈ pushforward_density(full, exp, 1.3; h = 0.15) rtol = 1e-10
end
