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
  check_point_estimator((s, x) -> conditional_variance_standard_error(s, x; h = 0.2), a, hist; points = [-1.0, 0.0, 1.0], outside = 100.0)

  @testset "conditional_variance_standard_error — the residual approximation (§4.2)" begin
    hbw = 0.4
    # positivity, scalar == vector, exact n^{-1/2} copy-merge scaling
    check_standard_error((s, x) -> conditional_variance_standard_error(s, x; h = hbw), hist;
                         points = [-1.0, 0.0, 1.0])
    # paired reader matches the separate calls; the SE is common to both variance methods
    r = conditional_variance_with_standard_error(hist, 0.5; h = hbw)
    @test r.value == conditional_variance(hist, 0.5; h = hbw)
    @test r.standard_error == conditional_variance_standard_error(hist, 0.5; h = hbw)
    @test conditional_variance_with_standard_error(hist, 0.5; h = hbw, method = :plugin).standard_error ==
          conditional_variance_with_standard_error(hist, 0.5; h = hbw, method = :residual).standard_error
    @test conditional_variance_with_standard_error(hist, [-0.5, 0.5]; h = hbw).standard_error ==
          conditional_variance_standard_error(hist, [-0.5, 0.5]; h = hbw)
    # ... and each half against its closed form: σ²(x) = 1 − ρ², and the residual SE √(ℓ(x)(μ₄ − σ⁴))
    # with the Gaussian μ₄ = 3σ⁴ and the interior leverage ℓ(x) ≈ R(K)/(N h f_X(x)), R(K) = 1/(2√π)
    σ2 = 1 - ρ^2
    for x in (0.0, 0.5)
      fX = exp(-x^2 / 2) / sqrt(2π)
      closed = sqrt(2 * σ2^2 / (2 * sqrt(π) * N * hbw * fX))
      rx = conditional_variance_with_standard_error(hist, x; h = hbw)
      @test rx.value ≈ σ2 rtol = 0.05
      @test rx.standard_error ≈ closed rtol = 0.15
    end
    @test isnan(conditional_variance_standard_error(hist, 100.0; h = hbw))
    @test_throws ArgumentError conditional_variance_standard_error(hist, 0.0; h = -1.0)
    @test_throws ArgumentError conditional_variance_standard_error(hist, 0.0; h = hbw, degree = 3)

    # the residual approximation captures conditional KURTOSIS: a heavy-tailed (Student-t₅) response, at the
    # same σ²(x), has a larger μ₄ and so a wider variance band than the Gaussian oracle
    Random.seed!(84); pert = 300_000
    t5(n) = randn(n) ./ sqrt.((randn(n) .^ 2 .+ randn(n) .^ 2 .+ randn(n) .^ 2 .+ randn(n) .^ 2 .+ randn(n) .^ 2) ./ 5)
    Xg = randn(pert); heavy = Histogram2D(-6.0, 6.0, 160, -8.0, 8.0, 320)
    observe_pairs!(heavy, Xg, ρ .* Xg .+ sqrt(1 - ρ^2) .* (t5(pert) ./ sqrt(5 / 3)))
    Xn = randn(pert); light = Histogram2D(-6.0, 6.0, 160, -8.0, 8.0, 320)
    observe_pairs!(light, Xn, ρ .* Xn .+ sqrt(1 - ρ^2) .* randn(pert))
    @test conditional_variance_standard_error(heavy, 0.5; h = hbw) > conditional_variance_standard_error(light, 0.5; h = hbw)

    # Monte-Carlo: the SE predicts the empirical sampling spread of σ̂²(0.5) across independent datasets
    Random.seed!(85); K = 30; per = 120_000
    chunks = Histogram2D[]; ests = Float64[]
    for _ in 1:K
      Xk = randn(per); Yk = ρ .* Xk .+ sqrt(1 - ρ^2) .* randn(per)
      s = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120); observe_pairs!(s, Xk, Yk)
      push!(chunks, s); push!(ests, conditional_variance(s, 0.5; h = hbw))
    end
    @test conditional_variance_standard_error(chunks[1], 0.5; h = hbw) ≈ std(ests) rtol = 0.2

    # coverage: the 95% band σ̂²(x) ± 1.96·SE brackets the constant 1−ρ² across the interior (aggregate; the
    # residual SE is a variance band, so the check lives in the variance-dominated regime, moderate n — the
    # fraction runs ≈ 0.95 across seeds here, the linear-mean oracle carrying negligible variance bias)
    xcov = -1.25:0.25:1.25
    covered = count(x -> abs(conditional_variance(chunks[1], x; h = hbw) - (1 - ρ^2)) ≤
                         1.96 * conditional_variance_standard_error(chunks[1], x; h = hbw), xcov)
    @test covered / length(xcov) ≥ 0.7
  end
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
  # the order-2 Y-power-sums, so n_y is taken fine.
  Random.seed!(31)
  N = 1_500_000
  ρ = 0.6
  v = 1 - ρ^2
  X = randn(N)
  Y = ρ .* X .+ sqrt(v) .* randn(N)
  g = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  observe_pairs!(g, X, Y)
  h_x = 0.4

  @testset "⟨g(Y)|x⟩ oracle — E[cos Y|x], E[sin Y|x]" begin
    for x in (-1.0, 0.0, 1.0)
      @test conditional_expectation(g, cos, x; h = h_x) ≈ cos(ρ * x) * exp(-v / 2) atol = 0.025
    end
    @test conditional_expectation(g, sin, 0.8; h = h_x) ≈ sin(ρ * 0.8) * exp(-v / 2) atol = 0.025
    @test conditional_expectation(g, sin, 0.0; h = h_x) ≈ 0.0 atol = 0.01
  end

  @testset "conditional_mean is exactly conditional_expectation(identity)" begin
    for x in (-1.5, -0.5, 0.5, 1.5)
      @test conditional_mean(g, x; h = h_x) == conditional_expectation(g, identity, x; h = h_x)
      @test conditional_mean(g, x; h = h_x, degree = 2) ==
            conditional_expectation(g, identity, x; h = h_x, degree = 2)
    end
  end

  @testset "plug-in variance = ⟨Y²|x⟩ − ⟨Y|x⟩² via the general reader" begin
    for x in (-1.0, 0.0, 1.0)
      plug = conditional_variance(g, x; h = h_x, method = :plugin)
      viaexp = conditional_expectation(g, y -> y^2, x; h = h_x) -
               conditional_expectation(g, identity, x; h = h_x)^2
      @test plug ≈ viaexp rtol = 1e-10
    end
  end

  @testset "scalar == vector; argument errors" begin
    xs = [-1.0, 0.0, 1.0]
    @test conditional_expectation(g, cos, xs; h = h_x) ≈
          [conditional_expectation(g, cos, x; h = h_x) for x in xs]
    @test_throws ArgumentError conditional_expectation(g, cos, 0.0; h = -1.0)
    @test_throws ArgumentError conditional_expectation(g, cos, 0.0; h = h_x, degree = 3)
  end

  # the universal query-point contract (scalar==vector, functional of the merged summary, NaN off support)
  half = N ÷ 2
  a = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  b = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 500)
  observe_pairs!(a, view(X, 1:half), view(Y, 1:half))
  observe_pairs!(b, view(X, half+1:N), view(Y, half+1:N))
  merge!(a, b)
  check_point_estimator((s, x) -> conditional_expectation(s, cos, x; h = h_x), a, g;
                        points = [-1.0, 0.0, 1.0], outside = 100.0)
  check_point_estimator((s, x) -> conditional_expectation_standard_error(s, cos, x; h = h_x), a, g;
                        points = [-1.0, 0.0, 1.0], outside = 100.0)
  check_standard_error((s, x) -> conditional_expectation_standard_error(s, cos, x; h = h_x), g;
                       points = [-1.0, 0.0, 1.0])
  check_se_bandwidth_scaling((s, x, hh) -> conditional_expectation_standard_error(s, identity, x; h = hh),
                             g, 0.0; h = 0.1, rtol = 0.06)

  @testset "conditional_pushforward: Z = exp(Y)|x ~ Lognormal(ρx, 1−ρ²)" begin
    x0 = 0.5
    lognormal_pdf(z) = 1 / (z * sqrt(2π * v)) * exp(-(log(z) - ρ * x0)^2 / (2v))
    for z in (0.8, 1.3, 2.0)
      @test conditional_pushforward_density(g, exp, x0, z; h_x = h_x, h_z = 0.18) ≈
            lognormal_pdf(z) rtol = 0.08
    end
    # exact normalization ∫ f̂ dz = 1 (symmetric kernel_z)
    zs = range(0.02, 15.0; length = 1800)
    curve = conditional_pushforward_density(g, exp, x0, collect(zs); h_x = h_x, h_z = 0.18)
    integral = sum(curve) * step(zs) - 0.5 * (curve[1] + curve[end]) * step(zs)
    @test integral ≈ 1.0 atol = 0.04
    # scalar == vector, and the h_z > 0 guard
    @test conditional_pushforward_density(g, exp, x0, [0.8, 1.3]; h_x = h_x, h_z = 0.18) ≈
          [conditional_pushforward_density(g, exp, x0, z; h_x = h_x, h_z = 0.18) for z in (0.8, 1.3)]
    @test_throws ArgumentError conditional_pushforward_density(g, exp, x0, 1.0; h_x = h_x, h_z = -0.1)
    @test_throws ArgumentError conditional_pushforward_density(g, exp, x0, 1.0; h_x = -1.0, h_z = 0.2)

    # standard error (leverage × the response-kernel variance): paired consistency, MC, off-support, throws
    rp = conditional_pushforward_density_with_standard_error(g, exp, x0, 1.0; h_x = h_x, h_z = 0.18)
    @test rp.value == conditional_pushforward_density(g, exp, x0, 1.0; h_x = h_x, h_z = 0.18)
    @test rp.standard_error == conditional_pushforward_density_standard_error(g, exp, x0, 1.0; h_x = h_x, h_z = 0.18)
    @test conditional_pushforward_density_with_standard_error(g, exp, x0, [0.8, 1.3]; h_x = h_x, h_z = 0.18).standard_error ==
          conditional_pushforward_density_standard_error(g, exp, x0, [0.8, 1.3]; h_x = h_x, h_z = 0.18)
    # ... and each half against its closed form: the conditional lognormal pdf, and the two-term variance
    # √(R(K)/(N h_x f_X(x)) · (R(K) f/h_z − f²)) — the interior local-linear leverage R(K)/(N h_x f_X(x)) times
    # the response-kernel variance R(K) f/h_z − f², R(K) = 1/(2√π), f_X the N(0, 1) marginal
    fz = lognormal_pdf(1.0); fX0 = exp(-x0^2 / 2) / sqrt(2π); R_K = 1 / (2 * sqrt(π))
    @test rp.value ≈ fz rtol = 0.08
    @test rp.standard_error ≈ sqrt(R_K / (N * h_x * fX0) * (R_K * fz / 0.18 - fz^2)) rtol = 0.12
    @test isnan(conditional_pushforward_density_standard_error(g, exp, 100.0, 1.0; h_x = h_x, h_z = 0.18))
    @test_throws ArgumentError conditional_pushforward_density_standard_error(g, exp, x0, 1.0; h_x = h_x, h_z = -0.1)
    @test_throws ArgumentError conditional_pushforward_density_standard_error(g, exp, x0, 1.0; h_x = h_x, h_z = 0.18, degree = 3)
    # 200 datasets put the spread's own error near 5%; 12% separates the two-term variance from the form
    # with −f²/n in place of the leverage, which runs ~17% high here
    Random.seed!(83); chunks = Histogram2D[]; mc = Float64[]
    for _ in 1:200
      Xk = randn(80_000); Yk = ρ .* Xk .+ sqrt(1 - ρ^2) .* randn(80_000)
      s = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 200); observe_pairs!(s, Xk, Yk)
      push!(chunks, s); push!(mc, conditional_pushforward_density(s, exp, x0, 1.0; h_x = h_x, h_z = 0.18))
    end
    @test conditional_pushforward_density_standard_error(chunks[1], exp, x0, 1.0; h_x = h_x, h_z = 0.18) ≈ std(mc) rtol = 0.12
    # coverage: the 95% band brackets the truth f_{exp(Y)|x}(z) across the interior. As for a 2-D density the
    # SE is a variance band and the estimate's bias dominates it at large n, so coverage is checked in the
    # variance-dominated regime (small n, here 5000) and is an aggregate interior fraction, not a pointwise gate.
    Random.seed!(96)
    Xc = randn(5000); Yc = ρ .* Xc .+ sqrt(1 - ρ^2) .* randn(5000)
    sc = Histogram2D(-6.0, 6.0, 200, -6.0, 6.0, 200); observe_pairs!(sc, Xc, Yc)
    zcov = 0.6:0.2:2.4
    covered = count(z -> abs(conditional_pushforward_density(sc, exp, x0, z; h_x = h_x, h_z = 0.18) - lognormal_pdf(z)) ≤
                         1.96 * conditional_pushforward_density_standard_error(sc, exp, x0, z; h_x = h_x, h_z = 0.18), zcov)
    @test covered / length(zcov) ≥ 0.6
  end

  @testset "coincidence: pushforward(identity, degree = 0) ≈ conditional_density" begin
    # both estimate f_{Y|X}; they agree up to the within-Y-cell correction (pushforward uses the degree-2
    # y-quadrature, conditional_density the cell-center value), i.e. to O(Δ_y²).
    for y in (-0.5, 0.3, 1.0)
      pf = conditional_pushforward_density(g, identity, 0.5, y; h_x = h_x, h_z = 0.2, degree = 0)
      cd = conditional_density(g, 0.5, y; h_x = h_x, h_y = 0.2)
      @test pf ≈ cd rtol = 0.02
    end
    # the standard errors coincide too at degree 0 (the pushforward SE reduces to the conditional-density SE)
    for y in (-0.5, 0.3, 1.0)
      @test conditional_pushforward_density_standard_error(g, identity, 0.5, y; h_x = h_x, h_z = 0.2, degree = 0) ≈
            conditional_density_standard_error(g, 0.5, y; h_x = h_x, h_y = 0.2) rtol = 0.03
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

  @testset "coverage: the 95% band brackets m(x) = ρx across the interior" begin
    # complement to the calibration above: |m̂(x) − ρx| ≤ 1.96·SE at each interior x (aggregate; local-linear
    # is unbiased for the linear mean, so per-point coverage sits near the nominal 0.95 — the interior fraction
    # runs ≈ 0.9 across seeds, dipping only when a whole-fit fluctuation shifts several correlated points at
    # once, so the floor is a lenient 0.6). SE via the identity reader.
    xcov = -1.5:0.25:1.5
    covered = count(x -> abs(conditional_mean(chunk, x; h = h) - ρ * x) ≤
                         1.96 * conditional_expectation_standard_error(chunk, identity, x; h = h), xcov)
    @test covered / length(xcov) ≥ 0.6
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

  @testset "conditional_expectation_with_standard_error — paired value + error, one fit" begin
    for x in (-1.0, 0.0, 1.0)
      r = conditional_expectation_with_standard_error(chunk, cos, x; h = h)
      @test r.value == conditional_expectation(chunk, cos, x; h = h)
      @test r.standard_error == conditional_expectation_standard_error(chunk, cos, x; h = h)
    end
    rv = conditional_expectation_with_standard_error(chunk, identity, [-1.0, 0.0, 1.0]; h = h)
    @test rv.value == conditional_expectation(chunk, identity, [-1.0, 0.0, 1.0]; h = h)
    @test rv.standard_error == conditional_expectation_standard_error(chunk, identity, [-1.0, 0.0, 1.0]; h = h)
    @test isnan(conditional_expectation_with_standard_error(chunk, identity, 100.0; h = h).standard_error)
    # ... and each half against its closed form: Y | x ∼ N(ρx, 1 − ρ²) gives E[cos Y | x] = cos(ρx) e^{-v/2}
    # and Var(cos Y | x) = (1 + cos(2ρx) e^{-2v})/2 − cos²(ρx) e^{-v}, v = 1 − ρ², and the SE is
    # √(R(K) Var/(n h f_X(x))) at the interior leverage, R(K) = 1/(2√π). The value carries the local-linear
    # O(h²) bias ½h² m''(x) ≈ −0.02 at h = 0.4, hence its absolute tolerance.
    v = 1 - ρ^2
    for x in (0.0, 0.5)
      Varx = (1 + cos(2ρ * x) * exp(-2v)) / 2 - cos(ρ * x)^2 * exp(-v)
      fX = exp(-x^2 / 2) / sqrt(2π)
      rx = conditional_expectation_with_standard_error(chunk, cos, x; h = h)
      @test rx.value ≈ cos(ρ * x) * exp(-v / 2) atol = 0.03
      @test rx.standard_error ≈ sqrt(Varx / (2 * sqrt(π) * per * h * fX)) rtol = 0.15
    end
    @test_throws ArgumentError conditional_expectation_with_standard_error(chunk, identity, 0.0; h = -1.0)
  end
end
