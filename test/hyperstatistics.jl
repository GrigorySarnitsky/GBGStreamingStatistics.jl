using Test
using GBGStreamingStatistics
using Random
include("testoracles.jl"); using .TestOracles

# A small builder: a fresh Float-count Histogram1D over the same data (deterministic, so the cone laws
# can be checked by reconstructing each side rather than copying).
_mk(data; lower = -8.0, upper = 8.0, nbins = 800) =
  (h = Histogram1D(lower, upper, nbins; count = Float64); observe!(h, data); h)


@testset "cone: scale! (the scalar action, §10)" begin
  # Float counts are required — the default UInt64 cannot hold a fractional weighted count
  @test_throws ArgumentError scale!(Histogram1D(-8.0, 8.0, 800), 0.5)

  Random.seed!(101)
  s = _mk(randn(100_000))
  m0, v0u, cm4, n0 = mean(s), var(s; corrected = false), central_moment(s, 4), sum(s.counts)
  scale!(s, 3.0)
  # scaling every count by a constant leaves the *normalized* distribution unchanged: the mean and the
  # UNcorrected central moments are exact invariants; the total mass scales by the weight
  @test mean(s) ≈ m0
  @test var(s; corrected = false) ≈ v0u
  @test central_moment(s, 4) ≈ cm4
  @test sum(s.counts) ≈ 3 * n0
  # (the corrected variance is *not* scale-invariant — its Bessel n/(n-1) legitimately shifts with n)
  @test var(s; corrected = true) ≈ v0u * (3n0) / (3n0 - 1)

  # homogeneity: scale by a then b equals scale by a·b
  Random.seed!(102); d = randn(60_000)
  twostep = _mk(d); scale!(twostep, 2.0); scale!(twostep, 3.0)
  onestep = _mk(d); scale!(onestep, 6.0)
  @test twostep.counts ≈ onestep.counts
  @test twostep.powersums ≈ onestep.powersums

  # weight = 1 is the identity
  Random.seed!(103); d = randn(40_000)
  s1 = _mk(d); scale!(s1, 1)
  @test s1.counts == _mk(d).counts
  # (scale! / weighted merge! are in-place broadcasts — their allocation-free audit is in allocations.jl)
end


@testset "cone: weighted merge! (the addition, §10)" begin
  Random.seed!(111); da = randn(50_000); db = randn(40_000) .+ 1.0

  # weight = 1 reduces to the plain monoid merge!
  i_w1 = _mk(da); merge!(i_w1, _mk(db); weight = 1)
  i_plain = _mk(da); merge!(i_plain, _mk(db))
  @test i_w1.counts == i_plain.counts && i_w1.powersums == i_plain.powersums

  # a weighted merge needs a Float-count `into`
  @test_throws ArgumentError merge!(Histogram1D(-8.0, 8.0, 800), _mk(da); weight = 0.5)
  @test (merge!(Histogram1D(-8.0, 8.0, 800), Histogram1D(-8.0, 8.0, 800); weight = 1); true)  # weight=1 ok on UInt64

  # weighted merge equals scale-then-merge:  into += w·from
  w = 0.37
  lhs = _mk(da); merge!(lhs, _mk(db); weight = w)
  rhs = _mk(da); fb = _mk(db); scale!(fb, w); merge!(rhs, fb)
  @test lhs.counts ≈ rhs.counts && lhs.powersums ≈ rhs.powersums

  # distributivity of the scalar action over the addition:  scale(a,w) + scale(b,w) = scale(a+b, w)
  left = _mk(da); scale!(left, w); tb = _mk(db); scale!(tb, w); merge!(left, tb)
  right = _mk(da); merge!(right, _mk(db)); scale!(right, w)
  @test left.counts ≈ right.counts && left.powersums ≈ right.powersums
end


@testset "hyperstatistics: the parameter-mixture density P(x) = Σ wₖ fₖ (§10)" begin
  # a two-component Gaussian mixture ½N(−2,1) + ½N(2,1), assembled by weighted merge of per-component
  # summaries — the discrete (and exact) form of the superstatistical superposition
  Random.seed!(121)
  T1 = Histogram1D(-10.0, 10.0, 1000); observe!(T1, randn(200_000) .- 2.0)
  T2 = Histogram1D(-10.0, 10.0, 1000); observe!(T2, randn(200_000) .+ 2.0)
  mix(x; w1 = 0.5) = w1 * gaussian_pdf(x; mu = -2.0) + (1 - w1) * gaussian_pdf(x; mu = 2.0)

  M = Histogram1D(-10.0, 10.0, 1000; count = Float64)
  merge!(M, T1; weight = 0.5 / sum(T1.counts))     # wₖ/Nₖ — each component normalized to a probability measure
  merge!(M, T2; weight = 0.5 / sum(T2.counts))
  @test sum(M.counts) ≈ 1.0                         # total mass = Σ wₖ
  # the mixture's mass is ≈ 1, so the *default* bandwidth (read off n ≈ 1) would oversmooth — supply h
  for x in (-3.0, -2.0, 0.0, 2.0, 3.0)
    @test density(M, x; h = 0.12) ≈ mix(x)  atol = 1e-2
  end

  # POST-HOC reweighting: the SAME per-component summaries, a different mixing law, no re-accumulation
  M73 = Histogram1D(-10.0, 10.0, 1000; count = Float64)
  merge!(M73, T1; weight = 0.7 / sum(T1.counts))
  merge!(M73, T2; weight = 0.3 / sum(T2.counts))
  for x in (-2.0, 0.0, 2.0)
    @test density(M73, x; h = 0.12) ≈ mix(x; w1 = 0.7)  atol = 1e-2
  end

  # equal component sizes ⇒ the bare weights wₖ keep the total a real sample count, so the default
  # bandwidth works directly (Nₖ is then a global constant that cancels in the density normalization)
  Mbare = Histogram1D(-10.0, 10.0, 1000; count = Float64)
  merge!(Mbare, T1; weight = 0.5); merge!(Mbare, T2; weight = 0.5)
  @test sum(Mbare.counts) ≈ Float64(sum(T1.counts))     # a representative n, not 1
  @test density(Mbare, 0.0) ≈ mix(0.0)  atol = 1e-2
end


@testset "hyperstatistics: χ²-superstatistics → Student-t (Beck & Cohen 2003)" begin
  # the canonical superstatistics result: X|β ~ N(0, 1/β) with the inverse-variance β ~ Gamma(ν/2, 2/ν)
  # has the heavy-tailed Student-t_ν marginal. Pool Gaussians over a β-grid with the Gamma mixing law.
  ν = 5.0
  Random.seed!(131)
  βs = range(0.06, 6.0; length = 60); Δβ = step(βs)
  M = Histogram1D(-15.0, 15.0, 3000; count = Float64)
  for β in βs
    w = gamma_pdf(β; shape = ν / 2, scale = 2 / ν) * Δβ
    w < 1e-7 && continue
    Tk = Histogram1D(-15.0, 15.0, 3000); observe!(Tk, (1 / sqrt(β)) .* randn(30_000))
    merge!(M, Tk; weight = w / sum(Tk.counts))
  end
  @test sum(M.counts) ≈ 1.0  atol = 2e-3            # ∫ f(β) dβ = 1

  # recovers the Student-t marginal through the body and into the heavy tail
  for x in (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)
    @test density(M, x; h = 0.10) ≈ student_t_pdf(x; nu = ν)  rtol = 5e-2
  end
  # the superstatistical tail is far heavier than any single Gaussian — the signature of the mixture
  @test density(M, 4.0; h = 0.10) > 20 * gaussian_pdf(4.0)
end


@testset "exponential forgetting (scale! ∘ observe!, §10)" begin
  # scaling the summary by λ < 1 before folding in the next block down-weights the past geometrically;
  # the recovered mean is then the EXACT λ-weighted mean of the two blocks (the shift formula is exact)
  Random.seed!(141)
  dA = randn(80_000) .- 2.0                          # an "old" regime
  dB = randn(50_000) .+ 3.0                          # a "new" regime
  λ = 0.4
  s = Histogram1D(-12.0, 12.0, 1200; count = Float64)
  observe!(s, dA); scale!(s, λ); observe!(s, dB)
  analytic = (λ * sum(dA) + sum(dB)) / (λ * length(dA) + length(dB))
  @test mean(s) ≈ analytic  rtol = 1e-10
  @test sum(s.counts) ≈ λ * length(dA) + length(dB)
end


@testset "Histogram2D cone (scale!, weighted merge!, §10)" begin
  # the joint summary lives in the same cone; the conditional estimators, being ratios of linear
  # functionals, are invariant under a scale! and carry mixture weights automatically
  Random.seed!(151)
  ρ = 0.6
  x = randn(300_000); y = ρ .* x .+ sqrt(1 - ρ^2) .* randn(300_000)
  h2 = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; count = Float64)
  observe_pairs!(h2, x, y)
  m_before = conditional_mean(h2, 1.0; h = 0.3)
  scale!(h2, 5.0)
  @test conditional_mean(h2, 1.0; h = 0.3) ≈ m_before       # ratio of linear functionals — scale cancels
  @test sum(h2.counts) ≈ 5 * 300_000                         # mass scaled

  # the structural cone laws hold cell-wise, exactly as for Histogram1D
  @test_throws ArgumentError scale!(Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120), 0.5)
  a = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; count = Float64); observe_pairs!(a, x, y)
  w1 = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; count = Float64); observe_pairs!(w1, x, y); merge!(w1, a; weight = 1)
  plain = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120; count = Float64); observe_pairs!(plain, x, y); merge!(plain, a)
  @test w1.counts == plain.counts && w1.powersums == plain.powersums
end
