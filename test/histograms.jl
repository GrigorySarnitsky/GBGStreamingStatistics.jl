using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random, Statistics
const H = Histograms

@testset "Histogram1D" begin
  @testset "construction & binning" begin
    h = Histogram1D(-5.0, 5.0, 10)
    @test H._nbins(h) == 10
    @test H._order(h) == 4
    @test bin_index(h, -5.0) == 1
    @test bin_index(h, -4.9) == 1
    @test bin_index(h, 0.0) == 6
    @test bin_index(h, 4.99) == 10
    @test bin_index(h, 5.0) == 11     # upper is the open edge → out of range
    @test bin_index(h, -6.0) < 1      # underflow
    @test_throws ArgumentError Histogram1D(5.0, -5.0, 10)
    @test_throws ArgumentError Histogram1D(-5.0, 5.0, 0)
    @test_throws ArgumentError Histogram1D(-5.0, 5.0, 10; order = 1)
  end

  @testset "counts & out-of-range skipping" begin
    h = Histogram1D(0.0, 10.0, 10)
    observe!(h, [0.5, 0.6, 5.5, 9.9, -1.0, 100.0])   # last two are out of range → skipped
    @test H._total(h) == 4
    @test h.counts[1] == 2            # 0.5, 0.6
    @test h.counts[6] == 1            # 5.5
    @test h.counts[10] == 1           # 9.9
  end

  @testset "moment recovery — Gaussian oracle, exact vs direct" begin
    Random.seed!(20240624)
    data = randn(2_000_000)
    h = Histogram1D(-8.0, 8.0, 1600; order = 6)
    observe!(h, data)
    @test H._total(h) == count(x -> -8.0 ≤ x < 8.0, data)

    # Exact agreement with the direct sample central moments (the Pébay recursion is exact arithmetic).
    m = mean(data)
    direct(p) = sum((data .- m) .^ p) / length(data)
    for p in 2:6
      @test H._sample_central_moment(h, p) ≈ direct(p) rtol = 1e-8 atol = 1e-12
    end

    # Agreement with the true Gaussian moments, to sampling error.
    @test H._sample_mean(h) ≈ 0 atol = 5e-3
    @test H._sample_variance(h) ≈ 1 rtol = 1e-2
    @test abs(H._sample_central_moment(h, 3)) < 1e-2          # skewness ≈ 0
    @test H._sample_central_moment(h, 4) ≈ 3 rtol = 2e-2      # kurtosis
    @test abs(H._sample_central_moment(h, 5)) < 5e-2          # ≈ 0
    @test H._sample_central_moment(h, 6) ≈ 15 rtol = 5e-2
  end

  @testset "Float32 input accumulates in Float64" begin
    Random.seed!(3)
    h = Histogram1D(-8.0, 8.0, 1000; order = 4)
    observe!(h, randn(Float32, 300_000))
    @test eltype(h.powersums) === Float64
    @test eltype(h.counts) === UInt64
    @test H._sample_variance(h) ≈ 1 rtol = 3e-2
  end
end
