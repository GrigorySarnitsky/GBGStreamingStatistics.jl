using Test
using GBGStreamingStatistics
using Random
const H = GBGStreamingStatistics.Histograms
include("testoracles.jl"); using .TestOracles

@testset "Bandwidth (bandwidth selection)" begin
  # A known sample: σ̂ recovered exactly from the binned moments.
  h = Histogram1D(0.0, 10.0, 100)
  observe!(h, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
  n = H._total(h)
  σ = GBGStreamingStatistics.Expectations._sample_std(h)
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
  hln = Histogram1D(0.0, 12.0, 2400); observe!(hln, exp.(0.6 .* randn(150_000)))
  gln = collect(0.05:0.05:6.0)
  tln = lognormal_pdf.(gln; sigma = 0.6)
  iseln(b) = ise(density(hln, gln; kernel = gamma_kernel, h = b), tln, 0.05)
  @test iseln(bandwidth(hln; rule = :plugin, kernel = gamma_kernel)) <
        0.8 * iseln(bandwidth(hln; rule = :normal, kernel = gamma_kernel))

  # …and on gamma(5) data
  Random.seed!(413)
  hg5 = Histogram1D(0.0, 20.0, 2000); observe!(hg5, [sum(-log(rand()) for _ in 1:5) for _ in 1:150_000])
  gg5 = collect(0.1:0.1:16.0)
  tg5 = gamma_pdf.(gg5; shape = 5)
  iseg5(b) = ise(density(hg5, gg5; kernel = gamma_kernel, h = b), tg5, 0.1)
  @test iseg5(bandwidth(hg5; rule = :plugin, kernel = gamma_kernel)) <
        0.9 * iseg5(bandwidth(hg5; rule = :normal, kernel = gamma_kernel))   # seed-robust margin (audit)

  # the beta-kernel plug-in on Beta(2,2) data
  Random.seed!(414)
  betadata = [(s = -log(rand()) - log(rand()); t = -log(rand()) - log(rand()); s / (s + t)) for _ in 1:150_000]
  hb = Histogram1D(0.0, 1.0, 800); observe!(hb, betadata)
  gb = collect(0.02:0.02:0.98)
  tb = beta_pdf.(gb; a = 2, b = 2)
  iseb(b) = ise(density(hb, gb; kernel = beta_kernel, h = b), tb, 0.02)
  # (the smooth, near-symmetric Beta(2,2) is the modest case — the rule of thumb is already near-optimal
  # there, so across seeds the plug-in is *competitive* rather than decisively better, ≈ 0.85–0.94× by the
  # seed audit; the strong, seed-robust plug-in gains are the skewed gamma cases above. docs/TESTING.md)
  @test iseb(bandwidth(hb; rule = :plugin, kernel = beta_kernel)) <
        1.02 * iseb(bandwidth(hb; rule = :normal, kernel = beta_kernel))

  # near-exponential (sharp boundary peak): the hard case — the plug-in must not be materially worse than
  # the rule of thumb (it is ≈ neutral there)
  ge = collect(0.05:0.05:9.0)
  te = exponential_pdf.(ge)
  isee(b) = ise(density(he, ge; kernel = gamma_kernel, h = b), te, 0.05)
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
  @test_throws ArgumentError bandwidth(he; kernel = gamma_kernel, degree = 1)   # asymmetric takes no degree
end

@testset "Marron–Wand claw — multimodal stress (the field's standard hard case, §4.5)" begin
  # The claw (Marron & Wand 1992): a broad N(0,1) with five sharp σ = 0.1 spikes. The normal-reference
  # bandwidth (driven by the global σ̂) oversmooths the spikes into a single bump; the data-driven
  # selectors read the curvature off the bin counts and resolve them — the case ISJ/SJ exist for, and the
  # benchmark the field validates bandwidth selection on (here ISJ ≈ SJ; both crush the normal reference).
  Random.seed!(701)
  function claw(n)
    out = Vector{Float64}(undef, n)
    for i in 1:n
      u = rand()
      out[i] = u < 0.5 ? randn() : (floor(Int, (u - 0.5) * 10) / 2 - 1) + 0.1 * randn()
    end
    out
  end
  hc = Histogram1D(-6.0, 6.0, 3000); observe!(hc, claw(200_000))
  bn = bandwidth(hc; rule = :normal); bs = bandwidth(hc; rule = :sj); bi = bandwidth(hc; rule = :isj)
  @test bi < 0.4 * bn                          # ISJ resolves the spikes the normal reference oversmooths
  @test bs < 0.4 * bn                          # …and so does classic SJ
  g = collect(-3.0:0.01:3.0)
  iseclaw(bw) = ise(density(hc, g; h = bw), claw_pdf.(g), 0.01)
  @test iseclaw(bi) < 0.1 * iseclaw(bn)        # ≈ an order of magnitude lower ISE — the multimodal payoff
  @test iseclaw(bs) < 0.1 * iseclaw(bn)
end
