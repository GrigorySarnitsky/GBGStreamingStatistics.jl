using Test
using GBGStreamingStatistics
using Random, Statistics
include("testoracles.jl"); using .TestOracles

# A representative set of analytic PDFs (accuracy) and the n^{-2p/(2p+1)} convergence rate of the binned
# density across its knobs. The estimand is f_X (§4.1); each PDF is an independent analytic oracle
# (`TestOracles`).
@testset "density: analytic-PDF accuracy & convergence rate (§4.1)" begin
  @testset "representative analytic PDFs" begin
    # recover the pdf at interior points and integrate to ≈ 1, from one accumulation
    function recovers(sampler, pdf, lower, upper, nbins, pts, n, atol_pdf, atol_int; seed = 11, kw = (;))
      Random.seed!(seed)
      h = Histogram1D(lower, upper, nbins); observe!(h, sampler(n))
      for x in pts
        @test density(h, x; kw...) ≈ pdf(x)  atol = atol_pdf
      end
      grid = range(lower, upper; length = 2001)
      @test sum(density(h, collect(grid); kw...)) * step(grid) ≈ 1  atol = atol_int
    end

    # 1. Gaussian — smooth, light tail
    recovers(n -> randn(n), gaussian_pdf, -10.0, 10.0, 2400, -3.0:1.0:3.0, 400_000, 6e-3, 2e-2)

    # 2. Uniform(0,1) — bounded support; the interior recovers f ≡ 1 (the edges carry the KDE boundary bias)
    recovers(n -> rand(n), uniform_pdf, -1.0, 2.0, 900, 0.2:0.2:0.8, 400_000, 1.5e-2, 2e-2)

    # 3. Exponential(1) — exponential tail; degree-1 local-log is exact (log f = −x is linear)
    recovers(n -> -log.(rand(n)), exponential_pdf, -1.0, 14.0, 1500, 0.5:0.5:4.0, 400_000, 1.2e-2, 2e-2)
    recovers(n -> -log.(rand(n)), exponential_pdf, -1.0, 14.0, 1500, 0.5:0.5:4.0, 400_000, 1.2e-2, 3.5e-2; kw = (degree = 1,))
    # (degree-1 is exact for the exponential — both degrees recover it above; the sharper degree-1 tail
    # win is the Laplace test in `test/densityestimation.jl`, where the feature is strong enough to be
    # robust at every point. For the milder exponential the finite-n advantage is marginal.)

    # 4. Student-t (ν = 3) — heavy power-law tail (the intermittency case)
    tsamp(n) = [randn() / sqrt((randn()^2 + randn()^2 + randn()^2) / 3) for _ in 1:n]
    tpdf(x) = student_t_pdf(x; nu = 3)
    recovers(tsamp, tpdf, -25.0, 25.0, 2500, -3.0:1.0:3.0, 400_000, 1.2e-2, 2e-2)
    # the variable bandwidth tracks the heavy tail better than the fixed one: smaller integrated error
    # over the tail x ∈ [4, 12] (≈ 3× here) — the intermittency payoff
    Random.seed!(13); ht = Histogram1D(-25.0, 25.0, 2500); observe!(ht, tsamp(400_000))
    tailpts = 4.0:0.5:12.0
    @test sum(abs(density(ht, x; adaptive = true) - tpdf(x)) for x in tailpts) <
          sum(abs(density(ht, x) - tpdf(x)) for x in tailpts)

    # 5. Log-normal (σ = 0.5) — skewed, one-sided
    recovers(n -> exp.(0.5 .* randn(n)), x -> lognormal_pdf(x; sigma = 0.5), -1.0, 10.0, 1100, 0.5:0.5:3.0, 400_000, 1.5e-2, 2e-2)

    # 6. Bimodal mixture ½N(−2.5,1) + ½N(2.5,1) — multimodality (both peaks and the valley)
    bipdf(x) = 0.5 * gaussian_pdf(x; mu = 2.5) + 0.5 * gaussian_pdf(x; mu = -2.5)
    Random.seed!(14); hb = Histogram1D(-8.0, 8.0, 1600)
    observe!(hb, vcat(2.5 .+ randn(200_000), -2.5 .+ randn(200_000)))
    @test density(hb, 2.5) ≈ bipdf(2.5)  atol = 1e-2        # peak
    @test density(hb, -2.5) ≈ bipdf(-2.5) atol = 1e-2        # peak
    @test density(hb, 0.0) ≈ bipdf(0.0)  rtol = 0.2         # valley (low density)
    grid = range(-8.0, 8.0; length = 2001)
    @test sum(density(hb, collect(grid))) * step(grid) ≈ 1  atol = 2e-2
  end

  @testset "convergence rate: ISE ∝ n^(−2p/(2p+1)) per degree" begin
    # accumulate Gaussian data at growing n; at the rate-optimal (default) bandwidth the integrated
    # squared error ISE = ∫(f̂ − φ)² should fall as n^(−2p/(2p+1)) — fit the log–log slope and check it
    # against the theoretical exponent. Bias order p = 2 (degree 0/1), p = 4 (degree 2).
    function ise_slope(; degree, adaptive = false, kernel = gaussian_kernel)
      ns = (30_000, 100_000, 300_000, 1_000_000)
      logn = Float64[]; logise = Float64[]
      for n in ns
        ises = Float64[]
        for s in 1:2
          Random.seed!(1000 + s)
          h = Histogram1D(-10.0, 10.0, 2400); observe!(h, randn(n))
          grid = range(-4.0, 4.0; length = 401)
          fh = density(h, collect(grid); degree = degree, adaptive = adaptive, kernel = kernel)
          push!(ises, ise(fh, gaussian_pdf.(grid), step(grid)))
        end
        push!(logn, log(n)); push!(logise, log(mean(ises)))
      end
      meann = mean(logn); meani = mean(logise)
      return sum((logn .- meann) .* (logise .- meani)) / sum((logn .- meann) .^ 2)
    end

    s0 = ise_slope(degree = 0, adaptive = false)
    s1 = ise_slope(degree = 1, adaptive = false)
    s2 = ise_slope(degree = 2, adaptive = false)
    sa = ise_slope(degree = 0, adaptive = true)
    se = ise_slope(degree = 0, kernel = epanechnikov_kernel)
    # the rate is kernel-independent (any second-order kernel ⇒ n^(−4/5)); the Epanechnikov, with its own
    # kernel-aware bandwidth, hits the same slope — a guard that the kernel-generic binned path delivers it
    # (only the constant differs, via δ(K), tested in test/bandwidth.jl). The non-Gaussian rate is a
    # deliberate omission: on a non-Gaussian the asymptotic regime is not reached by n = 10⁶ (the slope
    # muddies to ≈ −0.66 and degree 2 fails to separate from degree 0), so the Gaussian — whose fast
    # asymptotics make the rate visible at testable n — is the rate-test distribution.

    # degree 0/1 (bias O(h²)): MISE ∝ n^(−4/5) ⇒ slope ≈ −0.80 (RMISE ∝ n^(−2/5), the classic rate)
    @test s0 ≈ -0.80  atol = 0.10
    @test s1 ≈ -0.80  atol = 0.12
    # degree 2 (bias O(h⁴), the n^(−1/9) bandwidth): MISE ∝ n^(−8/9) ⇒ slope ≈ −0.89, and strictly
    # steeper than degree 0 — the higher-order rate the degree-aware bandwidth delivers
    @test s2 ≈ -0.89  atol = 0.10
    @test s2 < s0
    # the Abramson variable bandwidth shares the degree-0 order (it trades bias for variance pointwise,
    # not the asymptotic rate; whether the constant improves is density-dependent — it does not for this
    # smooth Gaussian, but does for peaked / heavy-tailed densities, per the ISE split commentary)
    @test sa ≈ -0.80  atol = 0.15
    @test se ≈ -0.80  atol = 0.12     # Epanechnikov — the same kernel-independent rate
  end

  @testset "binning-error convergence (WARPing → exact KDE): simple O(Δ), linear O(Δ²)" begin
    # The binned KDE places each bin's mass at a representative location; as the bins refine it converges
    # to the exact un-binned KDE (1/nh)Σ K((x−xᵢ)/h) at the same h. Simple binning (mass at the bin center)
    # is O(binwidth); linear binning (mass at the recovered within-bin centroid cⱼ + S₁[j]/nⱼ, which
    # matches the bin's first moment so the within-bin displacement averages to zero) is O(binwidth²) — the
    # default. Seed-averaged so the deterministic rate shows through the per-sample fluctuation.
    hbw = 0.3
    gridb = collect(-4.0:0.05:4.0)
    K(u) = exp(-u^2 / 2) / sqrt(2π)               # the smoothing kernel, as a reference
    nbs = (600, 1200, 2400, 4800)
    binwidths = 20 ./ nbs                          # range width 20 ([-10,10]) over nbins
    function discretization_errors(binning)
      errs = Float64[]
      for nb in nbs
        e = 0.0
        for sd in 1:3
          Random.seed!(sd); samples = randn(40_000)
          exact = [sum(K((x - s) / hbw) for s in samples) / (40_000 * hbw) for x in gridb]
          binned = Histogram1D(-10.0, 10.0, nb); observe!(binned, samples)
          e += sqrt(ise(density(binned, gridb; h = hbw, binning = binning), exact, 0.05))
        end
        push!(errs, e / 3)
      end
      errs
    end
    logslope(errs) = let lb = log.(binwidths), le = log.(errs), mb = mean(log.(binwidths)), me = mean(log.(errs))
      sum((lb .- mb) .* (le .- me)) / sum((lb .- mb) .^ 2)
    end
    errs_simple = discretization_errors(:simple)
    errs_linear = discretization_errors(:linear)
    @test errs_simple[end] < 0.4 * errs_simple[1]   # both converge as the grid refines
    @test errs_linear[end] < 0.4 * errs_linear[1]
    @test logslope(errs_simple) ≈ 1.0  atol = 0.4   # simple binning (mass at center) ⇒ O(binwidth)
    @test logslope(errs_linear) ≈ 2.0  atol = 0.4   # linear (centroid) binning ⇒ O(binwidth²)
    @test errs_linear[end] < 0.2 * errs_simple[end] # and far more accurate at a fine grid (the default)
  end
end
