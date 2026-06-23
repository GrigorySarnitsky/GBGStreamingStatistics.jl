using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random
const H = Histograms

@testset "mergeability (monoid laws)" begin
  Random.seed!(42)
  data = randn(600_000)
  mk() = Histogram1D(-8.0, 8.0, 800; order = 5)

  # accumulate-together
  full = mk(); observe!(full, data)

  # accumulate-then-merge over three independent chunks
  parts = [mk() for _ in 1:3]
  observe!(parts[1], view(data, 1:200_000))
  observe!(parts[2], view(data, 200_001:400_000))
  observe!(parts[3], view(data, 400_001:600_000))
  merged = mk()
  for p in parts
    merge!(merged, p)
  end

  @test merged.counts == full.counts                  # counts add — exact
  @test merged.means ≈ full.means rtol = 1e-9 atol = 1e-12
  for p in 2:5                                          # central moments — to round-off
    @test H._sample_central_moment(merged, p) ≈ H._sample_central_moment(full, p) rtol = 1e-8 atol = 1e-10
  end

  # merging an empty histogram is the identity
  before = copy(full.counts)
  before_var = H._sample_variance(full)
  merge!(full, mk())
  @test full.counts == before
  @test H._sample_variance(full) ≈ before_var

  # incompatible carriers are rejected
  @test_throws ArgumentError merge!(Histogram1D(-8.0, 8.0, 800), Histogram1D(-7.0, 8.0, 800))
  @test_throws ArgumentError merge!(Histogram1D(-8.0, 8.0, 800), Histogram1D(-8.0, 8.0, 400))
end
