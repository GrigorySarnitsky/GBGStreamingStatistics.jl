using Test
using GBGStreamingStatistics
using Random
using Statistics
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

  @testset "symmetric kernel options (parity): Epanechnikov in both directions" begin
    # The double-kernel takes a separate kernelx (conditioning) and kernely (response), each any ∫K=1
    # callable; Epanechnikov in both still integrates to 1 (the ∫K=1 normalization guarantees it) and
    # recovers the conditional to a few percent.
    ys = collect(-6.0:0.02:6.0)
    for x in (-0.5, 0.0, 0.5)
      fy = conditional_density(g, x, ys; hx = hx, hy = hy, kernelx = epanechnikov_kernel, kernely = epanechnikov_kernel)
      @test sum(fy) * 0.02 ≈ 1 atol = 1e-2
    end
    for x in (-0.5, 0.0, 0.5), y in (-0.5, 0.0, 0.5)
      @test conditional_density(g, x, y; hx = hx, hy = hy, kernelx = epanechnikov_kernel, kernely = epanechnikov_kernel) ≈
            cond_pdf(y, x) atol = 0.05
      @test transition_density(g, x, y; hx = hx, hy = hy, kernely = epanechnikov_kernel) ≈
            conditional_density(g, x, y; hx = hx, hy = hy, kernely = epanechnikov_kernel)
    end
  end

  @testset "asymmetric response kernel: bounded / positive Y (boundary-bias-free, §4.3)" begin
    # A bounded RESPONSE — Y ≥ 0 (a droplet radius). Independent oracle X ~ N(0,1), Y ~ Exponential(1),
    # so f(y|x) = e^{-y} for every x, positive at the boundary y = 0. A symmetric y-kernel leaks mass past
    # 0 and underestimates near the edge; the asymmetric gamma kernel (kernely = gamma_kernel) is
    # boundary-bias-free there. It is only ≈-normalized in y, so normalize = true makes ∫f̂dy = 1 exactly.
    Random.seed!(43)
    M = 2_000_000
    Xe = randn(M); Ye = -log.(rand(M))                          # Y ~ Exp(1) ⊥ X
    ge = Histogram2D(-6.0, 6.0, 240, 0.0, 12.0, 480; lag = 1)
    observe_pairs!(ge, Xe, Ye)
    fexp(y) = exp(-y)

    # boundary win: near y = 0 the gamma kernel beats the symmetric Gaussian, which leaks low
    for y0 in (0.1, 0.3)
      fg = conditional_density(ge, 0.0, y0; hx = 0.4, hy = 0.3, kernely = gamma_kernel, normalize = true)
      fs = conditional_density(ge, 0.0, y0; hx = 0.4, hy = 0.3)              # symmetric Gaussian
      @test abs(fg - fexp(y0)) < abs(fs - fexp(y0))           # gamma closer to the truth
    end
    @test conditional_density(ge, 0.0, 0.1; hx = 0.4, hy = 0.3) < fexp(0.1)  # Gaussian leaks low at the edge

    # interior recovery (the gamma kernel's width grows with y, so the convex e^{-y} is mildly oversmoothed
    # there — a looser tolerance than the boundary win), and exact normalization with the asymmetric kernel
    for y0 in (0.5, 1.0)
      @test conditional_density(ge, 0.0, y0; hx = 0.4, hy = 0.3, kernely = gamma_kernel, normalize = true) ≈
            fexp(y0) rtol = 0.12
    end
    ys = collect(0.0:0.02:12.0)
    @test sum(conditional_density(ge, 0.0, ys; hx = 0.4, hy = 0.3, kernely = gamma_kernel, normalize = true)) * 0.02 ≈
          1 atol = 2e-2
    # transition_density forwards kernely / normalize
    @test transition_density(ge, 0.0, 0.5; hx = 0.4, hy = 0.3, kernely = gamma_kernel, normalize = true) ≈
          conditional_density(ge, 0.0, 0.5; hx = 0.4, hy = 0.3, kernely = gamma_kernel, normalize = true)
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

  @testset "conditional_density_standard_error (double-kernel variability band, §4.3)" begin
    RK = 1 / (2 * sqrt(π))                                     # R(K) for the Gaussian
    fX(x) = gaussian_pdf(x; mu = 0.0, sigma = 1.0)            # the true marginal of X ∼ N(0,1)
    # oracle: the SE matches the two-term closed form √(R(K)²f(y|x)/(N hx hy f_X(x)) − f²/N) with the truth
    for (x0, y0) in ((0.0, 0.0), (0.5, 0.3))
      f = cond_pdf(y0, x0)
      closed = sqrt(RK^2 * f / (N * hx * hy * fX(x0)) - f^2 / N)
      @test conditional_density_standard_error(g, x0, y0; hx = hx, hy = hy) ≈ closed rtol = 0.06
    end
    # positivity, scalar == vector, and the exact n^{-1/2} copy-merge scaling
    check_standard_error((s, y) -> conditional_density_standard_error(s, 0.0, y; hx = hx, hy = hy), g;
                         points = [-0.3, 0.0, 0.3])
    # h^{-1/2} scaling in the response bandwidth (asymptotic; leading term at small hy, smooth point)
    check_se_bandwidth_scaling((s, y, hh) -> conditional_density_standard_error(s, 0.0, y; hx = hx, hy = hh),
                               g, 0.0; h = 0.15, rtol = 0.12)
    # paired reader: estimate and error match the separate calls exactly
    r = conditional_density_with_standard_error(g, 0.0, 0.2; hx = hx, hy = hy)
    @test r.value == conditional_density(g, 0.0, 0.2; hx = hx, hy = hy)
    @test r.standard_error == conditional_density_standard_error(g, 0.0, 0.2; hx = hx, hy = hy)
    rv = conditional_density_with_standard_error(g, 0.0, [-0.2, 0.0, 0.2]; hx = hx, hy = hy)
    @test rv.standard_error == conditional_density_standard_error(g, 0.0, [-0.2, 0.0, 0.2]; hx = hx, hy = hy)
    # transition-density aliases forward to the conditional readers
    @test transition_density_standard_error(g, 0.0, 0.0; hx = hx, hy = hy) ==
          conditional_density_standard_error(g, 0.0, 0.0; hx = hx, hy = hy)
    @test transition_density_with_standard_error(g, 0.0, 0.0; hx = hx, hy = hy).standard_error ==
          conditional_density_standard_error(g, 0.0, 0.0; hx = hx, hy = hy)
    # off-support NaN; bad bandwidth and the asymmetric-kernel variance (a recorded gap) throw
    @test isnan(conditional_density_standard_error(g, 100.0, 0.0; hx = hx, hy = hy))
    @test_throws ArgumentError conditional_density_standard_error(g, 0.0, 0.0; hx = -1.0, hy = hy)
    @test_throws ArgumentError conditional_density_standard_error(g, 0.0, 0.0; hx = hx, hy = hy, kernely = gamma_kernel)
    @test_throws ArgumentError conditional_density_standard_error(g, 0.0, 0.0; hx = hx, hy = hy, kernelx = gamma_kernel)

    # Monte-Carlo: the predicted SE matches the empirical spread of f̂(0|0) across independent datasets
    Random.seed!(71); K = 30; per = 80_000
    chunks = Histogram2D[]; ests = Float64[]
    for _ in 1:K
      Xk = randn(per); Yk = ρ .* Xk .+ sqrt(1 - ρ^2) .* randn(per)
      s = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1); observe_pairs!(s, Xk, Yk)
      push!(chunks, s); push!(ests, conditional_density(s, 0.0, 0.0; hx = hx, hy = hy))
    end
    @test conditional_density_standard_error(chunks[1], 0.0, 0.0; hx = hx, hy = hy) ≈ std(ests) rtol = 0.2

    # coverage: the 95% band f̂(y|x) ± 1.96·SE brackets the truth N(ρx, 1−ρ²) across an interior (x, y) grid.
    # A standard error is a *variance* band; the estimate carries an O(hx²+hy²) bias the band excludes, and for
    # a 2-D kernel density that bias dominates the (tight) SE at large n — so coverage is checked in the
    # variance-dominated regime (small n, here 5000) and is an aggregate interior fraction, never a strict
    # pointwise gate (the fraction runs ≈ 0.8 across seeds; the 0.6 floor only guards against a broken SE scale).
    Random.seed!(94)
    Xc = randn(5000); Yc = ρ .* Xc .+ sqrt(1 - ρ^2) .* randn(5000)
    sc = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1); observe_pairs!(sc, Xc, Yc)
    grid = [(x, y) for x in (-0.5, -0.25, 0.0, 0.25, 0.5), y in (-1.0, -0.5, 0.0, 0.5, 1.0)]
    covered = count(((x, y),) -> abs(conditional_density(sc, x, y; hx = hx, hy = hy) - cond_pdf(y, x)) ≤
                                1.96 * conditional_density_standard_error(sc, x, y; hx = hx, hy = hy), grid)
    @test covered / length(grid) ≥ 0.6
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
    check_point_estimator((s, y) -> conditional_density_standard_error(s, x0, y; hx = hx, hy = hy), ga, g; points = [-0.3, 0.0, 0.3], rtol = 1e-9)
  end
end
