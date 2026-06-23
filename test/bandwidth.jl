using Test
using GBGStreamingStatistics
const H = GBGStreamingStatistics.Histograms

@testset "Bandwidth (bandwidth selection)" begin
  # On a known sample, Silverman = (4/3)^(1/5)·σ̂·n^(-1/5), Scott = σ̂·n^(-1/5).
  h = Histogram1D(0.0, 10.0, 100)
  observe!(h, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
  n = H._total(h)
  σ = H._sample_std(h)
  @test σ ≈ sqrt(60 / 9)                      # exact recovery from the binned moments
  @test bandwidth(h) ≈ (4 / 3)^(1 / 5) * σ * float(n)^(-1 / 5)
  @test bandwidth(h; rule = :scott) ≈ σ * float(n)^(-1 / 5)
  @test bandwidth(h) > bandwidth(h; rule = :scott)

  # errors
  @test_throws ArgumentError bandwidth(h; rule = :nope)
  @test_throws ArgumentError bandwidth(Histogram1D(0.0, 1.0, 10))   # n = 0 < 2
  h1 = Histogram1D(0.0, 1.0, 10); observe!(h1, 0.5)
  @test_throws ArgumentError bandwidth(h1)                          # n = 1 < 2
end
