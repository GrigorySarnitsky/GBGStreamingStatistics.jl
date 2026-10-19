using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random, Statistics
const H = Histograms

# Recover the X/Y marginal moments and covariance from a Histogram2D's per-cell power-sums (the 2D
# shift formula) — exercises that the bivariate power-sums S_10 … S_02 are accumulated correctly.
function recover_2d(h)
  counts = h.counts; S = h.powersums; n_x = h.n_x
  n = 0.0; sx = 0.0; sy = 0.0
  for cell in 1:length(counts)
    nc = counts[cell]; iszero(nc) && continue
    i_x = (cell - 1) % n_x + 1; i_y = (cell - 1) ÷ n_x + 1
    n += nc
    sx += H._xbincenter(h, i_x) * nc + S[1, cell]
    sy += H._ybincenter(h, i_y) * nc + S[2, cell]
  end
  xbar = sx / n; ybar = sy / n
  M2x = 0.0; M2y = 0.0; Cxy = 0.0
  for cell in 1:length(counts)
    nc = counts[cell]; iszero(nc) && continue
    i_x = (cell - 1) % n_x + 1; i_y = (cell - 1) ÷ n_x + 1
    dx = H._xbincenter(h, i_x) - xbar; dy = H._ybincenter(h, i_y) - ybar
    M2x += nc * dx^2 + 2 * dx * S[1, cell] + S[3, cell]
    M2y += nc * dy^2 + 2 * dy * S[2, cell] + S[5, cell]
    Cxy += nc * dx * dy + dx * S[2, cell] + dy * S[1, cell] + S[4, cell]
  end
  return xbar, ybar, M2x / n, M2y / n, Cxy / n
end

@testset "Histogram2D" begin
  @testset "construction & binning" begin
    @test_throws ArgumentError Histogram2D(5.0, -5.0, 4, 0.0, 1.0, 4)
    @test_throws ArgumentError Histogram2D(0.0, 1.0, 0, 0.0, 1.0, 4)
    hd = Histogram2D(0.0, 4.0, 4, 0.0, 4.0, 4; lag = 3)
    @test hd.lag == 3
    observe_pairs!(hd, [0.5, 1.5, 3.5], [0.5, 0.5, 3.5])      # cells 1, 2, 16
    @test H._total(hd) == 3
    @test hd.counts[1] == 1 && hd.counts[2] == 1 && hd.counts[16] == 1
    @test_throws ArgumentError observe_pairs!(hd, [1.0], [1.0, 2.0])
  end

  @testset "independent-Gaussian oracle (marginals + covariance)" begin
    Random.seed!(202406)
    n = 1_000_000
    X = randn(n); Y = randn(n)                                # independent
    h = Histogram2D(-6.0, 6.0, 120, -6.0, 6.0, 120)
    observe_pairs!(h, X, Y)
    @test H._total(h) == count(i -> -6 ≤ X[i] < 6 && -6 ≤ Y[i] < 6, 1:n)
    xbar, ybar, vx, vy, cov = recover_2d(h)
    @test xbar ≈ 0 atol = 5e-3
    @test ybar ≈ 0 atol = 5e-3
    @test vx ≈ 1 rtol = 2e-2
    @test vy ≈ 1 rtol = 2e-2
    @test abs(cov) < 1e-2                                      # X ⟂ Y ⇒ covariance ≈ 0
  end

  @testset "mergeability + allocation" begin
    Random.seed!(11)
    X = randn(400_000); Y = randn(400_000)
    full = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(full, X, Y)
    summary_a = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(summary_a, view(X, 1:200_000), view(Y, 1:200_000))
    summary_b = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(summary_b, view(X, 200_001:400_000), view(Y, 200_001:400_000))
    merge!(summary_a, summary_b)
    @test summary_a.counts == full.counts
    @test summary_a.powersums ≈ full.powersums rtol = 1e-9 atol = 1e-10
    @test_throws ArgumentError merge!(Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80),
                                      Histogram2D(-5.0, 6.0, 80, -6.0, 6.0, 80))

    hb = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); xb = randn(10_000); yb = randn(10_000)
    observe_pairs!(hb, xb, yb)
    @test (@allocated observe_pairs!(hb, xb, yb)) < 4096       # bounded KA launch, no O(n) temporary
    ma = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(ma, xb, yb)
    mb = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(mb, xb, yb)
    merge!(ma, mb)
    @test (@allocated merge!(ma, mb)) == 0
  end
end
