using Test
using GBGStreamingStatistics
using Random
const H = GBGStreamingStatistics.Histograms

@testset "Bandwidth (bandwidth selection)" begin
  # A known sample: σ̂ recovered exactly from the binned moments.
  h = Histogram1D(0.0, 10.0, 100)
  observe!(h, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
  n = H._total(h)
  σ = H._sample_std(h)
  @test σ ≈ sqrt(60 / 9)                                # exact recovery from the binned moments

  # rule = :normal (default), Gaussian kernel, degree 0/1 → Silverman (4/3)^(1/5)·σ̂·n^(-1/5).
  @test bandwidth(h) ≈ (4 / 3)^(1 / 5) * σ * float(n)^(-1 / 5)  rtol = 1e-3
  @test bandwidth(h; degree = 1) ≈ bandwidth(h)  rtol = 1e-6           # deg 0,1 share the 2nd-order rule

  # the kernel enters only through the canonical factor δ(K): the Epanechnikov h is 2.214× the
  # Gaussian, i.e. ≈ 2.345·σ̂·n^(-1/5) (= statsmodels' normal_reference_constant).
  @test bandwidth(h; kernel = epanechnikov_kernel) ≈ 2.3449 * σ * float(n)^(-1 / 5)  rtol = 2e-3
  @test bandwidth(h; kernel = epanechnikov_kernel) / bandwidth(h) ≈ 2.2138  rtol = 2e-3

  # degree 2 is a higher-order (O(h⁴)) fit → the n^(-1/9) rate and a wider, more slowly shrinking h.
  # Matching the closed form tests both the constant and the -1/9 exponent.
  @test bandwidth(h; degree = 2) ≈ 1.0834 * σ * float(n)^(-1 / 9)  rtol = 3e-3
  @test bandwidth(h; degree = 2) > bandwidth(h; degree = 0)

  # errors
  @test_throws ArgumentError bandwidth(h; rule = :nope)
  @test_throws ArgumentError bandwidth(h; degree = -1)
  @test_throws ArgumentError bandwidth(h; degree = 5)                  # coverage gap: degree ≤ 3
  @test_throws ArgumentError bandwidth(Histogram1D(0.0, 1.0, 10))      # n = 0 < 2
  h1 = Histogram1D(0.0, 1.0, 10); observe!(h1, 0.5)
  @test_throws ArgumentError bandwidth(h1)                             # n = 1 < 2
end

@testset "Sheather–Jones plug-in (rule = :sj, §4.5)" begin
  # on Gaussian data the data-driven curvature ψ₄ ≈ the normal reference ⇒ SJ ≈ normal-reference
  Random.seed!(101)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(300_000))
  @test bandwidth(hg; rule = :sj) ≈ bandwidth(hg)  rtol = 5e-2

  # the eval kernel enters as the same canonical δ(K) rescale, exact (decoupled from the Gaussian pilot)
  @test bandwidth(hg; rule = :sj, kernel = epanechnikov_kernel) / bandwidth(hg; rule = :sj) ≈ 2.2138  rtol = 1e-2

  # on a bimodal mixture the normal reference oversmooths badly (σ̂ inflated by the inter-mode gap);
  # SJ sees the within-mode curvature and picks a much smaller h
  Random.seed!(102)
  hm = Histogram1D(-10.0, 10.0, 2000); observe!(hm, vcat(randn(150_000) .- 3, randn(150_000) .+ 3))
  @test bandwidth(hm; rule = :sj) < 0.5 * bandwidth(hm)

  # mergeable: SJ is a deterministic functional of the bin counts ⇒ accumulate-then-merge = together
  Random.seed!(103); d1 = randn(120_000); d2 = randn(120_000)
  a = Histogram1D(-8.0, 8.0, 1600); observe!(a, d1)
  b = Histogram1D(-8.0, 8.0, 1600); observe!(b, d2)
  both = Histogram1D(-8.0, 8.0, 1600); observe!(both, vcat(d1, d2))
  merge!(a, b)
  @test bandwidth(a; rule = :sj) ≈ bandwidth(both; rule = :sj)

  # Sheather–Jones is a second-order selector: degree ≥ 2 is a recorded coverage gap
  @test_throws ArgumentError bandwidth(hg; rule = :sj, degree = 2)
end

@testset "AIC fit-based selector (rule = :aic, §4.5)" begin
  # the AIC of the binned KDE has a stable interior optimum — mildly tighter than the normal reference
  # on Gaussian data (it adapts to the density rather than the worst-case reference)
  Random.seed!(201)
  hg = Histogram1D(-8.0, 8.0, 1600); observe!(hg, randn(250_000))
  @test 0.5 * bandwidth(hg) < bandwidth(hg; rule = :aic) < bandwidth(hg)

  # the kernel enters as the same exact canonical δ(K) rescale
  @test bandwidth(hg; rule = :aic, kernel = epanechnikov_kernel) / bandwidth(hg; rule = :aic) ≈ 2.2138  rtol = 1e-2

  # on a bimodal mixture AIC adapts to the within-mode curvature ⇒ far below the oversmoothing normal ref
  Random.seed!(202)
  hm = Histogram1D(-10.0, 10.0, 2000); observe!(hm, vcat(randn(125_000) .- 3, randn(125_000) .+ 3))
  @test bandwidth(hm; rule = :aic) < 0.5 * bandwidth(hm)

  # mergeable: AIC is a deterministic functional of the bin counts ⇒ accumulate-then-merge = together
  Random.seed!(203); d1 = randn(100_000); d2 = randn(100_000)
  a = Histogram1D(-8.0, 8.0, 1600); observe!(a, d1)
  b = Histogram1D(-8.0, 8.0, 1600); observe!(b, d2)
  both = Histogram1D(-8.0, 8.0, 1600); observe!(both, vcat(d1, d2))
  merge!(a, b)
  @test bandwidth(a; rule = :aic) ≈ bandwidth(both; rule = :aic)

  # degree-0/1 (KDE) selector: degree ≥ 2 is a coverage gap; the leave-one-out CV rules are not
  # provided at all (the binned-ties pathology) — they error as unknown rules
  @test_throws ArgumentError bandwidth(hg; rule = :aic, degree = 2)
  @test_throws ArgumentError bandwidth(hg; rule = :lscv)
  @test_throws ArgumentError bandwidth(hg; rule = :lcv)
end

@testset "Improved Sheather–Jones (rule = :isj, §4.5)" begin
  # on Gaussian data ISJ ≈ the normal reference (and ≈ classic SJ)
  Random.seed!(301)
  hg = Histogram1D(-8.0, 8.0, 2048); observe!(hg, randn(250_000))
  @test bandwidth(hg; rule = :isj) ≈ bandwidth(hg)  rtol = 5e-2
  # the eval kernel enters as the same exact canonical δ(K) rescale
  @test bandwidth(hg; rule = :isj, kernel = epanechnikov_kernel) / bandwidth(hg; rule = :isj) ≈ 2.2138  rtol = 1e-2

  # on a strongly multimodal mixture ISJ estimates the curvature nonparametrically (from the DCT) and
  # picks a much smaller h than the oversmoothing normal reference
  Random.seed!(302)
  hm = Histogram1D(-10.0, 10.0, 2048); observe!(hm, vcat(randn(125_000) .- 3, randn(125_000) .+ 3))
  @test bandwidth(hm; rule = :isj) < 0.5 * bandwidth(hm)

  # mergeable: ISJ is a deterministic functional of the bin counts (via their DCT)
  Random.seed!(303); d1 = randn(100_000); d2 = randn(100_000)
  a = Histogram1D(-8.0, 8.0, 2048); observe!(a, d1)
  b = Histogram1D(-8.0, 8.0, 2048); observe!(b, d2)
  both = Histogram1D(-8.0, 8.0, 2048); observe!(both, vcat(d1, d2))
  merge!(a, b)
  @test bandwidth(a; rule = :isj) ≈ bandwidth(both; rule = :isj)

  # degree ≥ 2 is a coverage gap
  @test_throws ArgumentError bandwidth(hg; rule = :isj, degree = 2)
end

@testset "Chen plug-in for asymmetric kernels (rule = :plugin, §4.5)" begin
  # rule = :normal is the σ̂ n^(−2/5) rule of thumb (the public Bessel-corrected `std` vs the internal
  # population σ̂ differ by ≈ √(n/(n−1)), hence rtol)
  Random.seed!(411)
  he = Histogram1D(0.0, 12.0, 2400); observe!(he, -log.(rand(150_000)))
  @test bandwidth(he; rule = :normal, kernel = gamma_kernel) ≈
        std(he) * float(sum(he.counts))^(-2 / 5)  rtol = 1e-3

  # the plug-in lowers the ISE on skewed positive data, where the rule of thumb mis-scales (the
  # droplet-radius regime). gamma kernel on log-normal data:
  Random.seed!(412)
  lnpdf(x) = exp(-(log(x))^2 / (2 * 0.36)) / (x * 0.6 * sqrt(2π))
  hln = Histogram1D(0.0, 12.0, 2400); observe!(hln, exp.(0.6 .* randn(150_000)))
  gln = collect(0.05:0.05:6.0)
  iseln(b) = sum(abs2, density(hln, gln; kernel = gamma_kernel, h = b) .- lnpdf.(gln)) * 0.05
  @test iseln(bandwidth(hln; rule = :plugin, kernel = gamma_kernel)) <
        0.8 * iseln(bandwidth(hln; rule = :normal, kernel = gamma_kernel))

  # …and on gamma(5) data
  Random.seed!(413)
  g5pdf(x) = x^4 * exp(-x) / 24
  hg5 = Histogram1D(0.0, 20.0, 2000); observe!(hg5, [sum(-log(rand()) for _ in 1:5) for _ in 1:150_000])
  gg5 = collect(0.1:0.1:16.0)
  iseg5(b) = sum(abs2, density(hg5, gg5; kernel = gamma_kernel, h = b) .- g5pdf.(gg5)) * 0.1
  @test iseg5(bandwidth(hg5; rule = :plugin, kernel = gamma_kernel)) <
        0.85 * iseg5(bandwidth(hg5; rule = :normal, kernel = gamma_kernel))

  # the beta-kernel plug-in on Beta(2,2) data
  Random.seed!(414)
  betadata = [(s = -log(rand()) - log(rand()); t = -log(rand()) - log(rand()); s / (s + t)) for _ in 1:150_000]
  hb = Histogram1D(0.0, 1.0, 800); observe!(hb, betadata)
  bpdf(x) = 6x * (1 - x); gb = collect(0.02:0.02:0.98)
  iseb(b) = sum(abs2, density(hb, gb; kernel = beta_kernel, h = b) .- bpdf.(gb)) * 0.02
  # (the smooth, near-symmetric Beta(2,2) is a modest-gain case — the strong gains are the skewed gamma
  # cases above; here the plug-in is reliably at least a little better)
  @test iseb(bandwidth(hb; rule = :plugin, kernel = beta_kernel)) <
        0.95 * iseb(bandwidth(hb; rule = :normal, kernel = beta_kernel))

  # near-exponential (sharp boundary peak): the hard case — the plug-in must not be materially worse than
  # the rule of thumb (it is ≈ neutral there)
  exppdf(x) = exp(-x); ge = collect(0.05:0.05:9.0)
  isee(b) = sum(abs2, density(he, ge; kernel = gamma_kernel, h = b) .- exppdf.(ge)) * 0.05
  @test isee(bandwidth(he; rule = :plugin, kernel = gamma_kernel)) <
        1.25 * isee(bandwidth(he; rule = :normal, kernel = gamma_kernel))
  @test bandwidth(he; rule = :plugin, kernel = gamma_kernel) > 0

  # mergeable: a deterministic functional of the bin counts
  Random.seed!(415); d1 = exp.(0.6 .* randn(80_000)); d2 = exp.(0.6 .* randn(80_000))
  a = Histogram1D(0.0, 12.0, 2400); observe!(a, d1)
  bh = Histogram1D(0.0, 12.0, 2400); observe!(bh, d2)
  both = Histogram1D(0.0, 12.0, 2400); observe!(both, vcat(d1, d2))
  merge!(a, bh)
  @test bandwidth(a; rule = :plugin, kernel = gamma_kernel) ≈
        bandwidth(both; rule = :plugin, kernel = gamma_kernel)  rtol = 1e-12

  # the symmetric selectors do not apply to an asymmetric kernel; degree ≥ 1 and :plugin-on-symmetric error
  @test_throws ArgumentError bandwidth(he; rule = :sj, kernel = gamma_kernel)
  @test_throws ArgumentError bandwidth(he; rule = :isj, kernel = gamma_kernel)
  @test_throws ArgumentError bandwidth(he; rule = :aic, kernel = gamma_kernel)
  @test_throws ArgumentError bandwidth(he; rule = :plugin, kernel = gamma_kernel, degree = 1)
  @test_throws ArgumentError bandwidth(he; rule = :plugin)   # :plugin needs an asymmetric kernel
end
