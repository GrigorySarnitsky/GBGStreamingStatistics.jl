using Test
using GBGStreamingStatistics
using Random

@testset "DensityEstimation (density)" begin
  K(u) = exp(-u^2 / 2) / sqrt(2π)

  # One point binned at a known bin center → the binned KDE is a single Gaussian bump at that center
  # (the WARPing approximation evaluates at bin centers, not the raw sample).
  h = Histogram1D(-1.0, 1.0, 2)        # bins [-1,0) center −0.5, [0,1) center 0.5
  observe!(h, 0.4)                      # lands in bin 2 (center 0.5)
  bw = 0.2
  @test density(h, 0.5; h = bw) ≈ K(0.0) / bw
  @test density(h, 0.3; h = bw) ≈ K((0.3 - 0.5) / bw) / bw

  # an empty histogram has zero density everywhere
  @test density(Histogram1D(-1.0, 1.0, 4), 0.0; h = 0.1) == 0.0

  # a vector of query points returns a vector
  ys = density(h, [0.3, 0.5, 0.7]; h = bw)
  @test ys isa AbstractVector
  @test length(ys) == 3
  @test ys[2] ≈ density(h, 0.5; h = bw)

  # the bandwidth must be positive
  @test_throws ArgumentError density(h, 0.0; h = 0.0)
  @test_throws ArgumentError density(h, 0.0; h = -0.5)
end

@testset "kernels: Gaussian + Epanechnikov" begin
  # the two shipped kernels are public
  @test gaussian_kernel(0.0) ≈ 1 / sqrt(2π)
  @test epanechnikov_kernel(0.0) == 0.75                  # 3/4 at the mode
  @test epanechnikov_kernel(0.5) == 0.75 * (1 - 0.25)
  @test epanechnikov_kernel(1.0) == 0.0                   # compact support [-1, 1]
  @test epanechnikov_kernel(2.0) == 0.0

  h = Histogram1D(-1.0, 1.0, 2); observe!(h, 0.4)         # one sample, bin 2 (center 0.5)
  bw = 0.6
  @test density(h, 0.5; h = bw, kernel = epanechnikov_kernel) ≈ epanechnikov_kernel(0.0) / bw
  @test density(h, 0.5; h = bw, kernel = gaussian_kernel) ≈ density(h, 0.5; h = bw)   # explicit == default
  @test density(h, -0.9; h = 0.3, kernel = epanechnikov_kernel) == 0.0                # compact support → 0 far away

  # both kernels integrate the recovered density to ≈ 1 over the range
  Random.seed!(7)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(200_000))
  xs = range(-8.0, 8.0; length = 3201); dx = step(xs)
  for kern in (gaussian_kernel, epanechnikov_kernel)
    @test sum(density(hg, x; h = 0.3, kernel = kern) for x in xs) * dx ≈ 1 atol = 2e-2
  end
end

@testset "density knobs: degree × fixed/variable bandwidth (§4.1)" begin
  Random.seed!(13)
  ha = Histogram1D(-8.0, 8.0, 1600); observe!(ha, randn(500_000))
  phi(x) = exp(-x^2 / 2) / sqrt(2π)

  # every corner of the grid recovers the Gaussian pdf, interior and into the tail
  for kw in ((;), (; adaptive = true), (; degree = 1), (; degree = 2), (; degree = 1, adaptive = true))
    @test density(ha, 0.0; kw...) ≈ phi(0.0) rtol = 5e-2
    @test density(ha, 3.0; kw...) ≈ phi(3.0) rtol = 1.5e-1
  end
  # degree-2 (log-quadratic) is exact for a Gaussian → tighter than the plain KDE
  @test density(ha, 1.0; degree = 2) ≈ phi(1.0) rtol = 2e-2

  # each variant integrates to ≈ 1 (local-likelihood need not be exactly 1)
  xs = range(-8.0, 8.0; length = 801); dx = step(xs)
  @test sum(density(ha, xs)) * dx ≈ 1 atol = 2e-2
  @test sum(density(ha, xs; adaptive = true)) * dx ≈ 1 atol = 2e-2
  @test sum(density(ha, xs; degree = 2)) * dx ≈ 1 atol = 2e-2
  @test sum(density(ha, xs; degree = 1)) * dx ≈ 1 atol = 3e-2

  # local-log-linear is exact for an exponential → it beats the KDE in a Laplace tail; and the variable
  # bandwidth, narrowing at the sharp Laplace cusp, resolves the peak the fixed KDE oversmooths. (The
  # adaptive payoff is for peaked / heavy-tailed densities — on a smooth Gaussian the wide tail kernel
  # instead oversmears the light tail and loses overall; see the ISE split in test/density_convergence.jl.)
  Random.seed!(2)
  hl = Histogram1D(-12.0, 12.0, 2400)
  observe!(hl, (-log.(rand(400_000))) .* rand((-1.0, 1.0), 400_000))   # Laplace(0, 1)
  flap(x) = 0.5 * exp(-abs(x))
  @test abs(density(hl, 5.0; degree = 1) - flap(5.0)) < abs(density(hl, 5.0) - flap(5.0))
  @test abs(density(hl, 0.0; adaptive = true) - flap(0.0)) < abs(density(hl, 0.0) - flap(0.0))   # cusp

  # h defaults to Silverman; vector form; empty; bad inputs
  @test density(ha, 0.0) ≈ density(ha, 0.0; h = bandwidth(ha))               # default h == Silverman
  ys = density(ha, [0.0, 1.0]; degree = 1, adaptive = true, kernel = epanechnikov_kernel)
  @test ys isa AbstractVector && length(ys) == 2 && all(>(0), ys)
  @test density(Histogram1D(-1.0, 1.0, 4), 0.0) == 0.0                       # empty → 0
  @test_throws ArgumentError density(ha, 0.0; h = -1.0)                      # bad bandwidth
  @test_throws ArgumentError density(ha, 0.0; degree = -1)                   # bad degree
end

@testset "density_standard_error (pointwise variability band, §4.1)" begin
  Random.seed!(31)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(400_000))
  n = 400_000; phi(x) = exp(-x^2 / 2) / sqrt(2π)
  RK = 1 / (2 * sqrt(π))                                  # R(K) for the Gaussian

  # the SE is exactly sqrt(R(K)·f̂/(n h))
  bw = bandwidth(hg)
  for x in (0.0, 1.0, 2.0)
    @test density_standard_error(hg, x; h = bw) ≈ sqrt(RK * density(hg, x; h = bw) / (n * bw))  rtol = 1e-6
  end

  # the 95% band f̂ ± 1.96·SE covers the truth across the interior (aggregate; bias dips it at mode/tail)
  xs = -2.5:0.25:2.5
  covered = count(x -> abs(density(hg, x; h = bw) - phi(x)) ≤ 1.96 * density_standard_error(hg, x; h = bw), xs)
  @test covered / length(xs) ≥ 0.75

  # more data ⇒ smaller SE
  Random.seed!(32); hs = Histogram1D(-8.0, 8.0, 1600); observe!(hs, randn(40_000))
  @test density_standard_error(hg, 0.0) < density_standard_error(hs, 0.0)

  # vector form; empty; bad inputs
  ses = density_standard_error(hg, [0.0, 1.0, 2.0])
  @test ses isa AbstractVector && length(ses) == 3 && all(>(0), ses)
  @test density_standard_error(Histogram1D(-1.0, 1.0, 4), 0.0) == 0.0
  @test_throws ArgumentError density_standard_error(hg, 0.0; h = -1.0)
  @test_throws ArgumentError density_standard_error(hg, 0.0; degree = 2)
end

@testset "boundary correction (bounds, §4.1)" begin
  # exponential on [0, ∞) (like droplet radii): the uncorrected KDE drops to ≈ half the true value at the
  # edge; the renormalized estimate (bounds = (0, Inf)) recovers it (Jones 1993)
  Random.seed!(41)
  h = Histogram1D(0.0, 12.0, 1200); observe!(h, -log.(rand(400_000)))
  expf(x) = exp(-x)
  @test abs(density(h, 0.02; bounds = (0.0, Inf)) - expf(0.02)) < abs(density(h, 0.02) - expf(0.02))
  @test density(h, 0.02; bounds = (0.0, Inf)) ≈ expf(0.02)  atol = 7e-2     # corrected ≈ true at the edge
  @test density(h, 0.02) < 0.7                                              # uncorrected ≈ half (true 0.98)

  # interior is unchanged (the kernel mass is ≈ 1 far from the edge)
  @test density(h, 2.0; bounds = (0.0, Inf)) ≈ density(h, 2.0)  rtol = 1e-6
  # the correction preserves normalization (integrates to ≈ 1 over the support)
  xs = range(0.0, 12.0; length = 2001)
  @test sum(density(h, collect(xs); bounds = (0.0, Inf))) * step(xs) ≈ 1  atol = 2e-2
  # kernel-general: Epanechnikov also lifts the boundary estimate
  @test density(h, 0.05; bounds = (0.0, Inf), kernel = epanechnikov_kernel) >
        density(h, 0.05; kernel = epanechnikov_kernel)

  # vector form; degree ≥ 1 + bounds is rejected (the local-log is already boundary-adaptive)
  @test density(h, [0.05, 2.0]; bounds = (0.0, Inf)) isa AbstractVector
  @test_throws ArgumentError density(h, 0.1; degree = 1, bounds = (0.0, Inf))
end

@testset "asymmetric kernels: gamma [0,∞), beta [0,1] (§4.1)" begin
  # gamma kernel on exponential(1) — boundary-free on the half-line (the droplet-radii case)
  Random.seed!(41)
  he = Histogram1D(0.0, 12.0, 1200); observe!(he, -log.(rand(400_000)))
  expf(x) = exp(-x)
  # recovers the density at the boundary, where the symmetric KDE collapses to ≈ half — and is closer
  @test density(he, 0.02; kernel = gamma_kernel) ≈ expf(0.02)  atol = 6e-2
  @test abs(density(he, 0.02; kernel = gamma_kernel) - expf(0.02)) < abs(density(he, 0.02) - expf(0.02))
  # interior agrees with the truth
  for x in (0.5, 1.0, 2.0)
    @test density(he, x; kernel = gamma_kernel) ≈ expf(x)  rtol = 5e-2
  end
  # integrates to ≈ 1
  xs = range(0.0, 12.0; length = 2001)
  @test sum(density(he, collect(xs); kernel = gamma_kernel)) * step(xs) ≈ 1  atol = 2e-2

  # beta kernel on Uniform(0,1) — boundary-free at both edges (true f ≡ 1)
  Random.seed!(42)
  hu = Histogram1D(0.0, 1.0, 500); observe!(hu, rand(400_000))
  for x in (0.02, 0.2, 0.5, 0.8, 0.98)
    @test density(hu, x; kernel = beta_kernel) ≈ 1.0  atol = 8e-2
  end

  # mergeable: a deterministic functional of the bin counts
  Random.seed!(43); d1 = -log.(rand(120_000)); d2 = -log.(rand(120_000))
  a = Histogram1D(0.0, 12.0, 1200); observe!(a, d1)
  b = Histogram1D(0.0, 12.0, 1200); observe!(b, d2)
  both = Histogram1D(0.0, 12.0, 1200); observe!(both, vcat(d1, d2))
  merge!(a, b)
  @test density(a, 0.5; kernel = gamma_kernel) ≈ density(both, 0.5; kernel = gamma_kernel)

  # the asymmetric kernel is its own estimator: degree / adaptive / bounds / bad-h are rejected
  @test_throws ArgumentError density(he, 0.5; kernel = gamma_kernel, degree = 1)
  @test_throws ArgumentError density(he, 0.5; kernel = gamma_kernel, adaptive = true)
  @test_throws ArgumentError density(he, 0.5; kernel = gamma_kernel, bounds = (0.0, Inf))
  @test_throws ArgumentError density(he, 0.5; kernel = gamma_kernel, h = -1.0)
end
