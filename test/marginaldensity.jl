using Test
using GBGStreamingStatistics
using Random
using SpecialFunctions: loggamma
include("testoracles.jl"); using .TestOracles
include("testcontracts.jl"); using .TestContracts

@testset "MarginalDensity (density)" begin
  K(u) = exp(-u^2 / 2) / sqrt(2π)

  # One sample → a single Gaussian bump. With linear (centroid) binning — the default — the bump sits at
  # the sample's exact location (the within-bin centroid of a lone sample IS the sample); with simple
  # binning it sits at the bin center (the WARPing approximation evaluates at bin centers).
  h = Histogram1D(-1.0, 1.0, 2)        # bins [-1,0) center −0.5, [0,1) center 0.5
  observe!(h, 0.4)                      # lands in bin 2 (center 0.5, centroid 0.4)
  bw = 0.2
  @test density(h, 0.4; h = bw) ≈ K(0.0) / bw                            # linear: bump at the sample 0.4
  @test density(h, 0.2; h = bw) ≈ K((0.2 - 0.4) / bw) / bw
  @test density(h, 0.5; h = bw, binning = :simple) ≈ K(0.0) / bw         # simple: bump at the center 0.5
  @test density(h, 0.3; h = bw, binning = :simple) ≈ K((0.3 - 0.5) / bw) / bw

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

@testset "density: universal contract (scalar==vector, mergeable)" begin
  # the universal query-point contract applied by construction to the marginal density — see
  # test/testcontracts.jl. No off-support sentinel: a density is defined everywhere (≈ 0 in the tails).
  Random.seed!(61)
  data = randn(300_000)
  full = Histogram1D(-8.0, 8.0, 1600); observe!(full, data)
  a = Histogram1D(-8.0, 8.0, 1600); observe!(a, view(data, 1:150_000))
  b = Histogram1D(-8.0, 8.0, 1600); observe!(b, view(data, 150_001:300_000))
  merge!(a, b)
  bw = bandwidth(full)
  check_point_estimator((s, x) -> density(s, x; h = bw),                 a, full; points = [-1.0, 0.0, 1.0])
  check_point_estimator((s, x) -> density(s, x; h = bw, degree = 1),     a, full; points = [-1.0, 0.0, 1.0])
  check_point_estimator((s, x) -> density_standard_error(s, x; h = bw),  a, full; points = [-1.0, 0.0, 1.0])
end

@testset "density: affine equivariance — f̂_Y(ax+b) = f̂_X(x)/a (§4.1)" begin
  # An exact, oracle-free structural invariant of the binned KDE: for Y = aX + b with the histogram range
  # and the bandwidth scaled the same way, the bin assignments are identical and the 1/(nh) rescales
  # cleanly, so f̂_Y(a x + b) = f̂_X(x)/a to floating point. It catches a stray normalization / σ-scaling
  # factor that a single-scale (σ = 1) oracle would miss, and guarantees unit-consistency — the estimate
  # of a droplet radius cannot depend on whether it is measured in μm or m.
  Random.seed!(91)
  X = randn(200_000)
  a = 7.3; b = -2.1
  hX = Histogram1D(-8.0, 8.0, 1600); observe!(hX, X)
  hY = Histogram1D(a * (-8.0) + b, a * 8.0 + b, 1600); observe!(hY, a .* X .+ b)   # range scaled + shifted alike
  hbw = 0.3
  for degree in (0, 1), x in (-1.5, 0.0, 1.0, 2.5)
    @test density(hY, a * x + b; h = a * hbw, degree = degree) ≈
          density(hX, x; h = hbw, degree = degree) / a  rtol = 1e-10
  end
end

@testset "kernels: Gaussian + Epanechnikov" begin
  # the two shipped kernels are public
  @test gaussian_kernel(0.0) ≈ 1 / sqrt(2π)
  @test epanechnikov_kernel(0.0) == 0.75                  # 3/4 at the mode
  @test epanechnikov_kernel(0.5) == 0.75 * (1 - 0.25)
  @test epanechnikov_kernel(1.0) == 0.0                   # compact support [-1, 1]
  @test epanechnikov_kernel(2.0) == 0.0

  h = Histogram1D(-1.0, 1.0, 2); observe!(h, 0.4)         # one sample, bin 2 (center 0.5, centroid 0.4)
  bw = 0.6
  @test density(h, 0.4; h = bw, kernel = epanechnikov_kernel) ≈ epanechnikov_kernel(0.0) / bw   # bump at the centroid 0.4
  @test density(h, 0.4; h = bw, kernel = gaussian_kernel) ≈ density(h, 0.4; h = bw)   # explicit == default
  @test density(h, -0.9; h = 0.3, kernel = epanechnikov_kernel) == 0.0                # compact support → 0 far away

  # both kernels integrate the recovered density to ≈ 1 over the range
  Random.seed!(7)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(200_000))
  xs = range(-8.0, 8.0; length = 3201); dx = step(xs)
  for kern in (gaussian_kernel, epanechnikov_kernel)
    @test sum(density(hg, x; h = 0.3, kernel = kern) for x in xs) * dx ≈ 1 atol = 2e-2
  end
end

@testset "density options: degree × fixed/variable bandwidth (§4.1)" begin
  Random.seed!(13)
  ha = Histogram1D(-8.0, 8.0, 1600); observe!(ha, randn(500_000))

  # every corner of the grid recovers the Gaussian pdf, interior and into the tail
  for kw in ((;), (; adaptive = true), (; degree = 1), (; degree = 2), (; degree = 1, adaptive = true))
    @test density(ha, 0.0; kw...) ≈ gaussian_pdf(0.0) rtol = 5e-2
    @test density(ha, 3.0; kw...) ≈ gaussian_pdf(3.0) rtol = 1.5e-1
  end
  # degree-2 (log-quadratic) is exact for a Gaussian → tighter than the plain KDE
  @test density(ha, 1.0; degree = 2) ≈ gaussian_pdf(1.0) rtol = 2e-2

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
  n = 400_000
  RK = 1 / (2 * sqrt(π))                                  # R(K) for the Gaussian

  # the SE is exactly the two-term KDE variance √(R(K)·f̂/(n h) − f̂²/n) (Wand & Jones)
  bw = bandwidth(hg)
  for x in (0.0, 1.0, 2.0)
    f = density(hg, x; h = bw)
    @test density_standard_error(hg, x; h = bw) ≈ sqrt(RK * f / (n * bw) - f^2 / n)  rtol = 1e-6
  end

  # the 95% band f̂ ± 1.96·SE covers the truth across the interior (aggregate; bias dips it at mode/tail)
  xs = -2.5:0.25:2.5
  covered = count(x -> abs(density(hg, x; h = bw) - gaussian_pdf(x)) ≤ 1.96 * density_standard_error(hg, x; h = bw), xs)
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


@testset "density_standard_error tracks the estimator's options + density_with_standard_error (§4.1)" begin
  Random.seed!(41)
  hg = Histogram1D(-8.0, 8.0, 800); observe!(hg, randn(200_000))
  bw = 0.3

  # adaptive: the wider local (balloon) bandwidth in the tail lowers the variance vs fixed-h
  @test density_standard_error(hg, 2.0; h = bw, adaptive = true) < density_standard_error(hg, 2.0; h = bw)
  @test isfinite(density_standard_error(hg, 0.0; h = bw, adaptive = true))

  # start (Hjort–Glad): leading variance unchanged, finite and positive
  @test density_standard_error(hg, 0.5; h = bw, start = normal_start) > 0

  # the asymmetric-kernel variance is a recorded coverage gap
  @test_throws ArgumentError density_standard_error(hg, 0.5; h = bw, kernel = gamma_kernel)

  # bounds (Jones renormalization) inflates the variance at the edge vs the un-renormalized KDE there
  he = Histogram1D(0.0, 12.0, 800); observe!(he, -log.(rand(200_000)))
  @test density_standard_error(he, 0.1; h = 0.25, bounds = (0.0, Inf)) >
        density_standard_error(he, 0.1; h = 0.25)

  # density_with_standard_error: one option spec ⇒ estimate and error match the separate calls exactly
  for x in (-1.0, 0.0, 1.3)
    r = density_with_standard_error(hg, x; h = bw, degree = 1, adaptive = true)
    @test r.value == density(hg, x; h = bw, degree = 1, adaptive = true)
    @test r.standard_error == density_standard_error(hg, x; h = bw, degree = 1, adaptive = true)
  end
  rv = density_with_standard_error(hg, [-1.0, 0.0, 1.0]; h = bw)
  @test rv.value == density(hg, [-1.0, 0.0, 1.0]; h = bw)
  @test rv.standard_error == density_standard_error(hg, [-1.0, 0.0, 1.0]; h = bw)

  # Monte-Carlo: the SE predicts the actual sampling spread of the estimate at the density peak
  Random.seed!(51); K = 25; per = 40_000
  summaries = Histogram1D[]; ests = Float64[]
  for _ in 1:K
    s = Histogram1D(-8.0, 8.0, 800); observe!(s, randn(per))
    push!(summaries, s); push!(ests, density(s, 0.0; h = bw))
  end
  @test density_standard_error(summaries[1], 0.0; h = bw) ≈ std(ests) rtol = 0.2
end

@testset "boundary correction (bounds, §4.1)" begin
  # exponential on [0, ∞) (like droplet radii): the uncorrected KDE drops to ≈ half the true value at the
  # edge; the renormalized estimate (bounds = (0, Inf)) recovers it (Jones 1993)
  Random.seed!(41)
  h = Histogram1D(0.0, 12.0, 1200); observe!(h, -log.(rand(400_000)))
  @test abs(density(h, 0.02; bounds = (0.0, Inf)) - exponential_pdf(0.02)) < abs(density(h, 0.02) - exponential_pdf(0.02))
  @test density(h, 0.02; bounds = (0.0, Inf)) ≈ exponential_pdf(0.02)  atol = 7e-2     # corrected ≈ true at the edge
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
  # recovers the density at the boundary, where the symmetric KDE collapses to ≈ half — and is closer
  @test density(he, 0.02; kernel = gamma_kernel) ≈ exponential_pdf(0.02)  atol = 6e-2
  @test abs(density(he, 0.02; kernel = gamma_kernel) - exponential_pdf(0.02)) < abs(density(he, 0.02) - exponential_pdf(0.02))
  # interior agrees with the truth
  for x in (0.5, 1.0, 2.0)
    @test density(he, x; kernel = gamma_kernel) ≈ exponential_pdf(x)  rtol = 5e-2
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

@testset "parametric start (Hjort–Glad) + normalization (§4.1)" begin
  # a CORRECT start removes the KDE's smoothing bias (the correction f/p ≡ 1): at a deliberately wide
  # bandwidth the plain KDE flattens the Gaussian peak (≈ 10% low) while the normal start recovers it
  Random.seed!(51)
  hg = Histogram1D(-8.0, 8.0, 800); observe!(hg, randn(50_000))
  grid = collect(-4.0:0.05:4.0)
  truth = gaussian_pdf.(grid)
  @test density(hg, 0.0; h = 0.5) < 0.95 * gaussian_pdf(0.0)                                # the KDE is biased
  @test density(hg, 0.0; h = 0.5, start = normal_start) ≈ gaussian_pdf(0.0)  rtol = 2e-2    # the start fixes it
  @test ise(density(hg, grid; h = 0.5, start = normal_start), truth, 0.05) <
        0.05 * ise(density(hg, grid; h = 0.5), truth, 0.05)

  # a MATCHING start on skewed positive data (the droplet-radius case): the log-normal start tracks the
  # log-normal far better than the symmetric KDE, decisively so in the right tail
  Random.seed!(52)
  σln = 0.6
  hln = Histogram1D(0.0, 25.0, 1250); observe!(hln, exp.(σln .* randn(200_000)))
  gln = collect(0.1:0.05:6.0)
  tln = lognormal_pdf.(gln; sigma = σln)
  @test ise(density(hln, gln; start = lognormal_start), tln, 0.05) <
        0.3 * ise(density(hln, gln), tln, 0.05)
  for x in (2.0, 3.0, 4.0)
    @test density(hln, x; start = lognormal_start) ≈ lognormal_pdf(x; sigma = σln)  rtol = 0.1
  end

  # exponential_start on Exponential(1): the one-parameter start (rate λ̂ = 1/μ̂, from the mean alone)
  # removes the smoothing bias the symmetric KDE carries at a wide bandwidth (it rounds off the convex
  # e^{−x}), and recovers the density in the interior
  Random.seed!(54)
  hexp = Histogram1D(-1.0, 16.0, 1700); observe!(hexp, -log.(rand(200_000)))
  gexp = collect(0.5:0.05:8.0)
  texp = exponential_pdf.(gexp)
  @test ise(density(hexp, gexp; h = 0.5, start = exponential_start), texp, 0.05) <
        0.55 * ise(density(hexp, gexp; h = 0.5), texp, 0.05)              # the start removes the wide-h bias
  for x in (1.0, 2.0, 3.0)
    @test density(hexp, x; h = 0.5, start = exponential_start) ≈ exponential_pdf(x)  rtol = 0.1
  end

  # inverse_gaussian_start on IG(μ, λ): the moment-fitted start (mean μ̂, shape λ̂ = μ̂³/σ̂²) recovers a
  # positive density with a sharp peak and an exponentially-decaying right tail; sampled via the
  # Michael–Schucany–Haas (1976) transformation
  rand_ig(μ, λ) = let y = randn()^2
    x = μ + μ^2 * y / (2λ) - (μ / (2λ)) * sqrt(4μ * λ * y + μ^2 * y^2)
    rand() ≤ μ / (μ + x) ? x : μ^2 / x
  end
  Random.seed!(55)
  μig, λig = 1.0, 3.0
  hig = Histogram1D(-0.5, 10.0, 1050); observe!(hig, [rand_ig(μig, λig) for _ in 1:200_000])
  gig = collect(0.2:0.05:5.0)
  tig = inverse_gaussian_pdf.(gig; mu = μig, lambda = λig)
  @test ise(density(hig, gig; start = inverse_gaussian_start), tig, 0.05) <
        ise(density(hig, gig), tig, 0.05)
  for x in (0.6, 1.0, 1.5, 2.0)
    @test density(hig, x; start = inverse_gaussian_start) ≈
          inverse_gaussian_pdf(x; mu = μig, lambda = λig)  rtol = 0.12
  end

  # `normalize` makes ∫f̂ = 1 EXACTLY over the support — a parametric start (and an asymmetric kernel) is
  # not self-normalizing
  fine = collect(range(0.0, 25.0; length = 5001)); Δ = 25.0 / 5000
  trapz(v) = (sum(v) - 0.5 * (v[1] + v[end])) * Δ
  @test abs(trapz(density(hln, fine; start = gamma_start)) - 1.0) > 1e-3            # the raw start is off
  @test trapz(density(hln, fine; start = gamma_start, normalize = true)) ≈ 1.0  atol = 1e-6
  @test trapz(density(hln, fine; kernel = gamma_kernel, normalize = true)) ≈ 1.0  atol = 1e-6

  # mergeable: the start parameters merge through the recovered moments, the reweighted counts add
  Random.seed!(53); d1 = exp.(σln .* randn(100_000)); d2 = exp.(σln .* randn(100_000))
  a = Histogram1D(0.0, 25.0, 1250); observe!(a, d1)
  b = Histogram1D(0.0, 25.0, 1250); observe!(b, d2)
  both = Histogram1D(0.0, 25.0, 1250); observe!(both, vcat(d1, d2))
  merge!(a, b)
  for x in (0.5, 1.5, 3.0)
    @test density(a, x; start = lognormal_start, normalize = true) ≈
          density(both, x; start = lognormal_start, normalize = true)  rtol = 1e-9
  end

  # a parametric start is the degree-0, fixed-h, symmetric-kernel estimator — the other options are rejected
  @test_throws ArgumentError density(hg, 0.0; start = normal_start, degree = 1)
  @test_throws ArgumentError density(hg, 0.0; start = normal_start, adaptive = true)
  @test_throws ArgumentError density(hg, 0.0; start = normal_start, bounds = (0.0, Inf))
  @test_throws ArgumentError density(hg, 0.0; start = normal_start, kernel = gamma_kernel)
  # a positive-support start needs the sample mean inside its support
  hneg = Histogram1D(-10.0, 2.0, 200); observe!(hneg, randn(2_000) .- 4.0)
  @test_throws ArgumentError density(hneg, 0.5; start = gamma_start)
end

@testset "custom kernels and starts: user-built, through density + bandwidth (§4.1)" begin
  # (1) a custom SYMMETRIC kernel is just a callable with ∫K = 1 — here the quartic (biweight)
  #     K(u) = (15/16)(1−u²)². The kernel-aware `bandwidth` reads its R(K), μ₂(K) off quadrature, so it
  #     flows through with no manual constants.
  Random.seed!(81)
  hq = Histogram1D(-8.0, 8.0, 800); observe!(hq, randn(200_000))
  quartic(u) = abs(u) < 1 ? (15 / 16) * (1 - u^2)^2 : 0.0
  @test bandwidth(hq; kernel = quartic) > 0                              # kernel-aware h, by quadrature
  @test density(hq, 0.0; kernel = quartic) ≈ gaussian_pdf(0.0)  rtol = 5e-2
  @test density(hq, 1.0; kernel = quartic) ≈ gaussian_pdf(1.0)  rtol = 8e-2

  # (2) a custom ASYMMETRIC kernel via the public `AsymmetricKernel` constructor: rebuild the gamma
  #     kernel from its contract and confirm it matches the built-in through both `density` and the
  #     `:plugin` bandwidth (same pdf + same AMISE weight ⇒ identical).
  Random.seed!(82)
  he = Histogram1D(0.0, 12.0, 1200); observe!(he, -log.(rand(300_000)))
  gamma_pdf_fn(c, x, b) =
    c > 0 ? (ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b; exp((ρ - 1) * log(c) - c / b - ρ * log(b) - loggamma(ρ))) : 0.0
  my_gamma = AsymmetricKernel(gamma_pdf_fn, (0.0, Inf), x -> x)
  for x in (0.1, 1.0, 3.0)
    @test density(he, x; kernel = my_gamma) ≈ density(he, x; kernel = gamma_kernel)  rtol = 1e-12
  end
  @test bandwidth(he; rule = :plugin, kernel = my_gamma) ≈ bandwidth(he; rule = :plugin, kernel = gamma_kernel)  rtol = 1e-12

  # the two-argument constructor omits the plug-in weight: `density` and the rule of thumb still work,
  # but `:plugin` fails loud rather than use a wrong w(x).
  my_gamma_noweight = AsymmetricKernel(gamma_pdf_fn, (0.0, Inf))
  @test density(he, 1.0; kernel = my_gamma_noweight) ≈ density(he, 1.0; kernel = gamma_kernel)  rtol = 1e-12
  @test bandwidth(he; rule = :normal, kernel = my_gamma_noweight) > 0
  @test_throws ArgumentError bandwidth(he; rule = :plugin, kernel = my_gamma_noweight)

  # (3) a custom `ParametricStart`: a Rayleigh start (one parameter σ, fit from the mean). On Rayleigh
  #     data the matching start beats the plain KDE, recovers the density, composes with `normalize`, and
  #     stays mergeable — the whole semiparametric path is open to a user model.
  rayleigh_start = ParametricStart(
    (mean, variance) -> (mean * sqrt(2 / π),),
    (x, p) -> x > 0 ? log(x) - 2 * log(p[1]) - x^2 / (2 * p[1]^2) : -Inf,
    (0.0, Inf))
  Random.seed!(83)
  σray = 1.5
  raydata = σray .* sqrt.(-2 .* log.(rand(200_000)))                    # inverse-CDF Rayleigh(σ) draws
  hr = Histogram1D(0.0, 12.0, 1200); observe!(hr, raydata)
  raypdf(x) = x > 0 ? (x / σray^2) * exp(-x^2 / (2σray^2)) : 0.0
  gray = collect(0.1:0.05:6.0)
  # at a deliberately wide bandwidth the KDE oversmooths the peak; the matching start removes that bias
  @test density(hr, σray; h = 0.7) < 0.9 * raypdf(σray)                               # the wide KDE is biased low
  @test density(hr, σray; h = 0.7, start = rayleigh_start) ≈ raypdf(σray)  rtol = 3e-2
  @test ise(density(hr, gray; h = 0.7, start = rayleigh_start), raypdf.(gray), 0.05) <
        0.5 * ise(density(hr, gray; h = 0.7), raypdf.(gray), 0.05)

  # mergeable: σ̂ merges through the recovered mean, the reweighted counts add
  Random.seed!(84)
  r1 = σray .* sqrt.(-2 .* log.(rand(100_000))); r2 = σray .* sqrt.(-2 .* log.(rand(100_000)))
  ra = Histogram1D(0.0, 12.0, 1200); observe!(ra, r1)
  rb = Histogram1D(0.0, 12.0, 1200); observe!(rb, r2)
  rboth = Histogram1D(0.0, 12.0, 1200); observe!(rboth, vcat(r1, r2))
  merge!(ra, rb)
  @test density(ra, 2.0; start = rayleigh_start, normalize = true) ≈
        density(rboth, 2.0; start = rayleigh_start, normalize = true)  rtol = 1e-9
end
