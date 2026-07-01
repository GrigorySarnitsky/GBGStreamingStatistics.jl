using Test
using GBGStreamingStatistics
using Random
include("testcontracts.jl"); using .TestContracts

@testset "ConditionalExpectations (conditional_mean / conditional_variance)" begin
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
    # The default residual-based (Fan–Yao) estimator fits E[(Y−m̂)²|x] = σ²(x) directly, removing the
    # curved-E[Y²|x] bias the plug-in carries; for this linear-mean oracle it recovers the constant.
    for x in (-1.0, 0.0, 1.0)
      @test conditional_variance(hist, x; h = 0.2) ≈ 1 - ρ^2 rtol = 0.1
    end
    @test conditional_variance(hist, 0.0; h = 0.2) > 0
  end

  @testset "degree and kernel options (parity with density)" begin   # the new local_fit-backed options
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


@testset "Fan–Yao conditional variance vs the plug-in (curved mean, §4.2)" begin
  # Heteroscedastic, curved-mean oracle: X ~ N(0,1), m(x) = x² (so m'' = 2), σ²(x) = 0.2 + 0.1x². The
  # plug-in Ê[Y²|x] − m̂² fits the quartic E[Y²|x] = x⁴ + σ²(x) with a local polynomial, so it is badly
  # biased away from the centre; the residual-based estimator fits E[(Y−m̂)²|x] = σ²(x) (far flatter).
  Random.seed!(2025)
  M = 3_000_000
  Xh = randn(M)
  σ2(x) = 0.2 + 0.1 * x^2
  Yh = Xh .^ 2 .+ sqrt.(σ2.(Xh)) .* randn(M)
  hh = Histogram2D(-6.0, 6.0, 240, -8.0, 45.0, 240)
  observe_pairs!(hh, Xh, Yh)

  @testset "degree-2 recovers the heteroscedastic σ²(x)" begin
    # degree 2 matches the quadratic mean (first-stage m̂ unbiased) and the quadratic σ²(x).
    for x in (-1.0, -0.5, 0.0, 0.5, 1.0)
      @test conditional_variance(hh, x; h = 0.3, degree = 2) ≈ σ2(x) rtol = 0.12
    end
  end

  @testset "residual-based beats the plug-in in the curved tail" begin
    for x in (-1.5, 1.5)
      fy = conditional_variance(hh, x; h = 0.3, method = :residual)   # the default
      pl = conditional_variance(hh, x; h = 0.3, method = :plugin)
      @test abs(fy - σ2(x)) < abs(pl - σ2(x))      # Fan–Yao closer to the truth
      @test pl > 1.5 * σ2(x)                        # the plug-in is badly inflated by the m'' curvature
    end
  end

  @test_throws ArgumentError conditional_variance(hh, 0.0; h = 0.3, method = :bogus)
end


@testset "asymmetric kernel at degree 0 — admitted by the unified weighting (§4.2)" begin
  # The shared kernel-weighting atom admits an asymmetric kernel at degree 0 (a Nadaraya–Watson ratio with a
  # boundary kernel) — what the old symmetric-only local_fit rejected with a MethodError. This is the
  # lens-consistency unlock (every degree-0 readout takes any kernel), NOT a recommended regression
  # estimator: the Chen gamma/beta kernels are a *density* boundary tool, and their asymmetry interacts with
  # the regression design density, so for a bounded conditioning variable degree-1 local-linear — which is
  # boundary-adaptive AND unbiased for a linear mean — is the tool. An asymmetric kernel has no jet, so
  # degree ≥ 1 with one is a clean ArgumentError.
  Random.seed!(77)
  M = 1_000_000
  Xb = -log.(rand(M)); Yb = Xb .+ 0.3 .* randn(M)              # X ~ Exp(1) ≥ 0, E[Y|X=x] = x
  gb = Histogram2D(0.0, 12.0, 240, -2.0, 14.0, 240); observe_pairs!(gb, Xb, Yb)

  # admitted and finite (was a MethodError), interior and at the bounded edge
  @test isfinite(conditional_mean(gb, 2.0; h = 0.4, degree = 0, kernel = gamma_kernel))
  @test isfinite(conditional_mean(gb, 0.2; h = 0.4, degree = 0, kernel = gamma_kernel))
  # degree-1 local-linear is the right tool for the bounded edge — unbiased for the linear mean
  @test conditional_mean(gb, 0.5; h = 0.4, degree = 1) ≈ 0.5 atol = 0.03

  # the one guard: an asymmetric kernel has no jet ⇒ degree ≥ 1 throws (cleanly, not a MethodError)
  @test_throws ArgumentError conditional_mean(gb, 1.0; h = 0.4, degree = 1, kernel = gamma_kernel)
  @test_throws ArgumentError conditional_variance(gb, 1.0; h = 0.4, degree = 2, kernel = gamma_kernel)
end


@testset "conditional_expectation & pushforward — the general ⟨g(Y)|x⟩ (§4.2)" begin
  # Bivariate-Gaussian: Y|X=x ~ N(ρx, 1−ρ²), so ⟨cos Y|x⟩ = cos(ρx)·e^{−(1−ρ²)/2} (and sin likewise) —
  # a clean analytic oracle for the general conditional expectation. The y-quadrature of a general g uses
  # the order-2 Y-power-sums, so ny is taken fine.
  Random.seed!(31)
  N = 1_500_000
  ρ = 0.6
  v = 1 - ρ^2
  X = randn(N)
  Y = ρ .* X .+ sqrt(v) .* randn(N)
  g = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  observe_pairs!(g, X, Y)
  hx = 0.4

  @testset "⟨g(Y)|x⟩ oracle — E[cos Y|x], E[sin Y|x]" begin
    for x in (-1.0, 0.0, 1.0)
      @test conditional_expectation(g, cos, x; h = hx) ≈ cos(ρ * x) * exp(-v / 2) atol = 0.025
    end
    @test conditional_expectation(g, sin, 0.8; h = hx) ≈ sin(ρ * 0.8) * exp(-v / 2) atol = 0.025
    @test conditional_expectation(g, sin, 0.0; h = hx) ≈ 0.0 atol = 0.01
  end

  @testset "conditional_mean is exactly conditional_expectation(identity)" begin
    for x in (-1.5, -0.5, 0.5, 1.5)
      @test conditional_mean(g, x; h = hx) == conditional_expectation(g, identity, x; h = hx)
      @test conditional_mean(g, x; h = hx, degree = 2) ==
            conditional_expectation(g, identity, x; h = hx, degree = 2)
    end
  end

  @testset "plug-in variance = ⟨Y²|x⟩ − ⟨Y|x⟩² via the general reader" begin
    for x in (-1.0, 0.0, 1.0)
      plug = conditional_variance(g, x; h = hx, method = :plugin)
      viaexp = conditional_expectation(g, y -> y^2, x; h = hx) -
               conditional_expectation(g, identity, x; h = hx)^2
      @test plug ≈ viaexp rtol = 1e-10
    end
  end

  @testset "scalar == vector; argument errors" begin
    xs = [-1.0, 0.0, 1.0]
    @test conditional_expectation(g, cos, xs; h = hx) ≈
          [conditional_expectation(g, cos, x; h = hx) for x in xs]
    @test_throws ArgumentError conditional_expectation(g, cos, 0.0; h = -1.0)
    @test_throws ArgumentError conditional_expectation(g, cos, 0.0; h = hx, degree = 3)
  end

  # the universal query-point contract (scalar==vector, functional of the merged summary, NaN off support)
  half = N ÷ 2
  a = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  b = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  observe_pairs!(a, view(X, 1:half), view(Y, 1:half))
  observe_pairs!(b, view(X, half+1:N), view(Y, half+1:N))
  merge!(a, b)
  check_point_estimator((s, x) -> conditional_expectation(s, cos, x; h = hx), a, g;
                        points = [-1.0, 0.0, 1.0], outside = 100.0)
  check_point_estimator((s, x) -> conditional_expectation_standard_error(s, cos, x; h = hx), a, g;
                        points = [-1.0, 0.0, 1.0], outside = 100.0)
  check_standard_error((s, x) -> conditional_expectation_standard_error(s, cos, x; h = hx), g;
                       points = [-1.0, 0.0, 1.0])

  @testset "conditional_pushforward: Z = exp(Y)|x ~ Lognormal(ρx, 1−ρ²)" begin
    x0 = 0.5
    lognormal_pdf(z) = 1 / (z * sqrt(2π * v)) * exp(-(log(z) - ρ * x0)^2 / (2v))
    for z in (0.8, 1.3, 2.0)
      @test conditional_pushforward_density(g, exp, x0, z; hx = hx, hz = 0.18) ≈
            lognormal_pdf(z) rtol = 0.08
    end
    # exact normalization ∫ f̂ dz = 1 (symmetric kernelz)
    zs = range(0.02, 15.0; length = 1800)
    curve = conditional_pushforward_density(g, exp, x0, collect(zs); hx = hx, hz = 0.18)
    integral = sum(curve) * step(zs) - 0.5 * (curve[1] + curve[end]) * step(zs)
    @test integral ≈ 1.0 atol = 0.04
    # scalar == vector, and the hz > 0 guard
    @test conditional_pushforward_density(g, exp, x0, [0.8, 1.3]; hx = hx, hz = 0.18) ≈
          [conditional_pushforward_density(g, exp, x0, z; hx = hx, hz = 0.18) for z in (0.8, 1.3)]
    @test_throws ArgumentError conditional_pushforward_density(g, exp, x0, 1.0; hx = hx, hz = -0.1)
  end

  @testset "coincidence: pushforward(identity, degree = 0) ≈ conditional_density" begin
    # both estimate f_{Y|X}; they agree up to the within-Y-cell correction (pushforward uses the degree-2
    # y-quadrature, conditional_density the cell-center value), i.e. to O(Δ_y²).
    for y in (-0.5, 0.3, 1.0)
      pf = conditional_pushforward_density(g, identity, 0.5, y; hx = hx, hz = 0.2, degree = 0)
      cd = conditional_density(g, 0.5, y; hx = hx, hy = 0.2)
      @test pf ≈ cd rtol = 0.02
    end
  end
end


@testset "conditional_expectation_standard_error — leverage × noise-scale (§4.2)" begin
  # The SE is √(leverage(x)·Var(g(Y)|x)) — the finite-sample sampling error of the local fit. The gold
  # standard: it predicts the *actual* spread of the estimate across independent datasets. Build K
  # independent chunks, and compare the predicted SE (from one chunk) to the empirical std across chunks.
  Random.seed!(808)
  ρ = 0.6
  K = 24
  per = 60_000
  h = 0.4
  summaries = Histogram2D[]
  ests0 = Float64[]
  for _ in 1:K
    X = randn(per)
    Y = ρ .* X .+ sqrt(1 - ρ^2) .* randn(per)
    s = Histogram2D(-6.0, 6.0, 160, -6.0, 6.0, 320)
    observe_pairs!(s, X, Y)
    push!(summaries, s)
    push!(ests0, conditional_mean(s, 0.0; h = h))
  end
  chunk = summaries[1]

  @testset "Monte-Carlo: predicted SE = empirical sampling spread (at the density peak)" begin
    @test conditional_expectation_standard_error(chunk, identity, 0.0; h = h) ≈ std(ests0) rtol = 0.2
    # leading-order closed form √(R(K)σ²/(n h f_X(0))) — R(gaussian) = 1/(2√π), f_X(0) = 1/√(2π)
    closed0 = sqrt((1 / (2 * sqrt(π))) * (1 - ρ^2) / (per * h * (1 / sqrt(2π))))
    @test conditional_expectation_standard_error(chunk, identity, 0.0; h = h) ≈ closed0 rtol = 0.15
  end

  @testset "properties — positive, tail growth, off-support, general g, scalar == vector" begin
    @test conditional_expectation_standard_error(chunk, identity, 2.0; h = h) >
          conditional_expectation_standard_error(chunk, identity, 0.0; h = h)   # sparser tail ⇒ larger SE
    @test isnan(conditional_expectation_standard_error(chunk, identity, 100.0; h = h))
    @test conditional_expectation_standard_error(chunk, cos, 0.5; h = h) > 0
    xs = [-0.5, 0.0, 0.5]
    @test conditional_expectation_standard_error(chunk, identity, xs; h = h) ≈
          [conditional_expectation_standard_error(chunk, identity, x; h = h) for x in xs]
    @test_throws ArgumentError conditional_expectation_standard_error(chunk, identity, 0.0; h = -1.0)
    @test_throws ArgumentError conditional_expectation_standard_error(chunk, identity, 0.0; h = h, degree = 3)
  end
end
