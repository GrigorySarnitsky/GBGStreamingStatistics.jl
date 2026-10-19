using Test
using GBGStreamingStatistics
using Random
using Statistics
include("testoracles.jl"); using .TestOracles
include("testcontracts.jl"); using .TestContracts

@testset "ConditionalDensity (f_{Y|X}, transition density §4.3)" begin
  # Bivariate-Gaussian fixture: Y = ρX + √(1−ρ²)·ε ⇒ Y | X=x ∼ 𝒩(ρx, 1−ρ²) — an exact analytic
  # conditional density to check normalization and pointwise recovery against.
  Random.seed!(11)
  n = 2_000_000
  ρ = 0.6
  X = randn(n)
  Y = ρ .* X .+ sqrt(1 - ρ^2) .* randn(n)
  joint = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  observe_pairs!(joint, X, Y)
  h_x = 0.3; h_y = 0.2

  cond_pdf(y, x) = gaussian_pdf(y; mu = ρ * x, sigma = sqrt(1 - ρ^2))   # 𝒩(ρx, 1−ρ²)

  @testset "normalization: ∫ f̂(y|x) dy ≈ 1" begin
    ys = collect(-6.0:0.02:6.0)
    for x in (-1.0, 0.0, 1.0)
      fy = conditional_density(joint, x, ys; h_x = h_x, h_y = h_y)
      @test sum(fy) * 0.02 ≈ 1 atol = 1e-2
    end
  end

  @testset "recovers 𝒩(ρx, 1−ρ²)" begin
    for x in (-0.5, 0.0, 0.5), y in (-0.5, 0.0, 0.5)
      @test conditional_density(joint, x, y; h_x = h_x, h_y = h_y) ≈ cond_pdf(y, x) atol = 0.03
    end
  end

  @testset "symmetric kernel options (parity): Epanechnikov in both directions" begin
    # The double-kernel takes a separate kernel_x (conditioning) and kernel_y (response), each any ∫K=1
    # callable; Epanechnikov in both still integrates to 1 (the ∫K=1 normalization guarantees it) and
    # recovers the conditional to a few percent.
    ys = collect(-6.0:0.02:6.0)
    for x in (-0.5, 0.0, 0.5)
      fy = conditional_density(joint, x, ys; h_x = h_x, h_y = h_y, kernel_x = epanechnikov_kernel, kernel_y = epanechnikov_kernel)
      @test sum(fy) * 0.02 ≈ 1 atol = 1e-2
    end
    for x in (-0.5, 0.0, 0.5), y in (-0.5, 0.0, 0.5)
      @test conditional_density(joint, x, y; h_x = h_x, h_y = h_y, kernel_x = epanechnikov_kernel, kernel_y = epanechnikov_kernel) ≈
            cond_pdf(y, x) atol = 0.05
      @test transition_density(joint, x, y; h_x = h_x, h_y = h_y, kernel_y = epanechnikov_kernel) ≈
            conditional_density(joint, x, y; h_x = h_x, h_y = h_y, kernel_y = epanechnikov_kernel)
    end
  end

  @testset "asymmetric response kernel: bounded / positive Y (boundary-bias-free, §4.3)" begin
    # A bounded RESPONSE — Y ≥ 0 (a droplet radius). Independent oracle X ~ 𝒩(0,1), Y ~ Exponential(1),
    # so f(y|x) = e^{-y} for every x, positive at the boundary y = 0. A symmetric y-kernel leaks mass past
    # 0 and underestimates near the edge; the asymmetric gamma kernel (kernel_y = gamma_kernel) is
    # boundary-bias-free there. It is only ≈-normalized in y, so normalize = true makes ∫f̂dy = 1 exactly.
    Random.seed!(43)
    n_e = 2_000_000
    Xe = randn(n_e); Ye = -log.(rand(n_e))                          # Y ~ Exp(1) ⊥ X
    ge = Histogram2D(-6.0, 6.0, 240, 0.0, 12.0, 480; lag = 1)
    observe_pairs!(ge, Xe, Ye)
    fexp(y) = exp(-y)

    # boundary win: near y = 0 the gamma kernel beats the symmetric Gaussian, which leaks low
    for y0 in (0.1, 0.3)
      fg = conditional_density(ge, 0.0, y0; h_x = 0.4, h_y = 0.3, kernel_y = gamma_kernel, normalize = true)
      fs = conditional_density(ge, 0.0, y0; h_x = 0.4, h_y = 0.3)              # symmetric Gaussian
      @test abs(fg - fexp(y0)) < abs(fs - fexp(y0))           # gamma closer to the truth
    end
    @test conditional_density(ge, 0.0, 0.1; h_x = 0.4, h_y = 0.3) < fexp(0.1)  # Gaussian leaks low at the edge

    # interior recovery (the gamma kernel's width grows with y, so the convex e^{-y} is mildly oversmoothed
    # there — a looser tolerance than the boundary win), and exact normalization with the asymmetric kernel
    for y0 in (0.5, 1.0)
      @test conditional_density(ge, 0.0, y0; h_x = 0.4, h_y = 0.3, kernel_y = gamma_kernel, normalize = true) ≈
            fexp(y0) rtol = 0.12
    end
    ys = collect(0.0:0.02:12.0)
    @test sum(conditional_density(ge, 0.0, ys; h_x = 0.4, h_y = 0.3, kernel_y = gamma_kernel, normalize = true)) * 0.02 ≈
          1 atol = 2e-2
    # transition_density forwards kernel_y / normalize
    @test transition_density(ge, 0.0, 0.5; h_x = 0.4, h_y = 0.3, kernel_y = gamma_kernel, normalize = true) ≈
          conditional_density(ge, 0.0, 0.5; h_x = 0.4, h_y = 0.3, kernel_y = gamma_kernel, normalize = true)
  end

  @testset "estimand-specific API: transition_density alias, NaN off the conditioning support, bandwidth errors" begin
    ys = [-0.5, 0.0, 0.5]
    fv = conditional_density(joint, 0.0, ys; h_x = h_x, h_y = h_y)
    @test transition_density(joint, 0.0, 0.0; h_x = h_x, h_y = h_y) ≈ fv[2]      # alias = conditional
    @test transition_density(joint, 0.0, ys; h_x = h_x, h_y = h_y) == fv
    @test isnan(conditional_density(joint, 100.0, 0.0; h_x = h_x, h_y = h_y))    # off-support conditioning x ⇒ NaN
    @test_throws ArgumentError conditional_density(joint, 0.0, 0.0; h_x = -1.0, h_y = h_y)
    @test_throws ArgumentError conditional_density(joint, 0.0, 0.0; h_x = h_x, h_y = 0.0)
  end

  @testset "conditional_density_standard_error (double-kernel variability band, §4.3)" begin
    R_K = 1 / (2 * sqrt(π))                                     # R(K) for the Gaussian
    fX(x) = gaussian_pdf(x; mu = 0.0, sigma = 1.0)            # the true marginal of X ∼ 𝒩(0,1)
    # oracle: the SE matches the two-term closed form √(R(K)/(N h_x f_X(x)) · (R(K) f(y|x)/h_y − f²)) with the
    # truth — the sum of the squared x-weights times the variance of one response-kernel term
    for (x0, y0) in ((0.0, 0.0), (0.5, 0.3))
      f = cond_pdf(y0, x0)
      closed = sqrt(R_K / (n * h_x * fX(x0)) * (R_K * f / h_y - f^2))
      @test conditional_density_standard_error(joint, x0, y0; h_x = h_x, h_y = h_y) ≈ closed rtol = 0.06
    end
    # positivity, scalar == vector, and the exact n^{-1/2} copy-merge scaling
    check_standard_error((summary, y) -> conditional_density_standard_error(summary, 0.0, y; h_x = h_x, h_y = h_y), joint;
                         points = [-0.3, 0.0, 0.3])
    # h^{-1/2} scaling in the response bandwidth (asymptotic): it holds where f·h_y/R(K) is small, so the −f²
    # term is negligible — in the tail at y = 1.5 (f ≈ 0.09), not at the peak, where that term is half the first
    check_se_bandwidth_scaling((summary, y, hh) -> conditional_density_standard_error(summary, 0.0, y; h_x = h_x, h_y = hh),
                               joint, 1.5; h = 0.15, rtol = 0.12)
    # paired reader: estimate and error match the separate calls exactly
    r = conditional_density_with_standard_error(joint, 0.0, 0.2; h_x = h_x, h_y = h_y)
    @test r.value == conditional_density(joint, 0.0, 0.2; h_x = h_x, h_y = h_y)
    @test r.standard_error == conditional_density_standard_error(joint, 0.0, 0.2; h_x = h_x, h_y = h_y)
    # ... and each half against its closed form: 𝒩(ρx, 1 − ρ²) and the two-term double-kernel variance
    for (x0, y0) in ((0.0, 0.0), (0.5, 0.3))
      f = cond_pdf(y0, x0)
      rx = conditional_density_with_standard_error(joint, x0, y0; h_x = h_x, h_y = h_y)
      @test rx.value ≈ f atol = 0.03
      @test rx.standard_error ≈ sqrt(R_K / (n * h_x * fX(x0)) * (R_K * f / h_y - f^2)) rtol = 0.06
    end
    rv = conditional_density_with_standard_error(joint, 0.0, [-0.2, 0.0, 0.2]; h_x = h_x, h_y = h_y)
    @test rv.standard_error == conditional_density_standard_error(joint, 0.0, [-0.2, 0.0, 0.2]; h_x = h_x, h_y = h_y)
    # transition-density aliases forward to the conditional readers
    @test transition_density_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y) ==
          conditional_density_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y)
    @test transition_density_with_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y).standard_error ==
          conditional_density_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y)
    # off-support NaN; bad bandwidth and the asymmetric-kernel variance (a recorded gap) throw
    @test isnan(conditional_density_standard_error(joint, 100.0, 0.0; h_x = h_x, h_y = h_y))
    @test_throws ArgumentError conditional_density_standard_error(joint, 0.0, 0.0; h_x = -1.0, h_y = h_y)
    @test_throws ArgumentError conditional_density_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y, kernel_y = gamma_kernel)
    @test_throws ArgumentError conditional_density_standard_error(joint, 0.0, 0.0; h_x = h_x, h_y = h_y, kernel_x = gamma_kernel)

    # Monte-Carlo: the predicted SE matches the empirical spread of f̂(0|0) across independent datasets. 200
    # datasets put the spread's own error near 5%, so 12% separates the two-term variance from the form with
    # −f²/n in place of the local weight, which runs ~25% high at this point.
    Random.seed!(71); K = 200; per = 80_000
    chunks = Histogram2D[]; ests = Float64[]
    for _ in 1:K
      Xk = randn(per); Yk = ρ .* Xk .+ sqrt(1 - ρ^2) .* randn(per)
      summary = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1); observe_pairs!(summary, Xk, Yk)
      push!(chunks, summary); push!(ests, conditional_density(summary, 0.0, 0.0; h_x = h_x, h_y = h_y))
    end
    @test conditional_density_standard_error(chunks[1], 0.0, 0.0; h_x = h_x, h_y = h_y) ≈ std(ests) rtol = 0.12

    # coverage: the 95% band f̂(y|x) ± 1.96·SE brackets the truth 𝒩(ρx, 1−ρ²) across an interior (x, y) grid.
    # A standard error is a *variance* band; the estimate carries an O(h_x²+h_y²) bias the band excludes, and for
    # a 2-D kernel density that bias dominates the (tight) SE at large n — so coverage is checked in the
    # variance-dominated regime (small n, here 5000) and is an aggregate interior fraction, never a strict
    # pointwise gate (the fraction runs ≈ 0.8 across seeds; the 0.6 floor only guards against a broken SE scale).
    Random.seed!(94)
    Xc = randn(5000); Yc = ρ .* Xc .+ sqrt(1 - ρ^2) .* randn(5000)
    sc = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1); observe_pairs!(sc, Xc, Yc)
    grid = [(x, y) for x in (-0.5, -0.25, 0.0, 0.25, 0.5), y in (-1.0, -0.5, 0.0, 0.5, 1.0)]
    covered = count(((x, y),) -> abs(conditional_density(sc, x, y; h_x = h_x, h_y = h_y) - cond_pdf(y, x)) ≤
                                1.96 * conditional_density_standard_error(sc, x, y; h_x = h_x, h_y = h_y), grid)
    @test covered / length(grid) ≥ 0.6
  end

  # the universal query-point contract (scalar==vector over y, a functional of the merged summary),
  # applied by construction across a range of conditioning x — see test/testcontracts.jl
  half = n ÷ 2
  ga = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  gb = Histogram2D(-6.0, 6.0, 240, -6.0, 6.0, 240; lag = 1)
  observe_pairs!(ga, view(X, 1:half), view(Y, 1:half))
  observe_pairs!(gb, view(X, half+1:n), view(Y, half+1:n))
  merge!(ga, gb)
  for x0 in (-0.5, 0.0, 0.5)
    check_point_estimator((summary, y) -> conditional_density(summary, x0, y; h_x = h_x, h_y = h_y), ga, joint; points = [-0.3, 0.0, 0.3], rtol = 1e-9)
    check_point_estimator((summary, y) -> conditional_density_standard_error(summary, x0, y; h_x = h_x, h_y = h_y), ga, joint; points = [-0.3, 0.0, 0.3], rtol = 1e-9)
  end
end
