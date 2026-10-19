using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
const H = Histograms

# A bare accumulator with no merge! method of its own — exercises the generic erroring fallback.
struct _DummyAccumulator <: GBGStreamingStatistics.AbstractAccumulator end

@testset "Accumulators interface + small exact checks" begin
  @test Histogram1D <: GBGStreamingStatistics.AbstractAccumulator
  @test Histogram2D <: GBGStreamingStatistics.AbstractAccumulator

  # Three points into a single bin — recovered mean/variance match the hand computation exactly.
  h = Histogram1D(0.0, 1.0, 1)
  observe!(h, [0.2, 0.4, 0.9])
  @test H._total(h) == 3
  μ = (0.2 + 0.4 + 0.9) / 3
  @test GBGStreamingStatistics.Expectations._sample_mean(h) ≈ μ
  @test GBGStreamingStatistics.Expectations._sample_variance(h) ≈ ((0.2 - μ)^2 + (0.4 - μ)^2 + (0.9 - μ)^2) / 3

  # merge! of two single-point histograms equals observing both points together (the n_B = 1 case).
  summary_a = Histogram1D(0.0, 1.0, 1); observe!(summary_a, 0.2)
  summary_b = Histogram1D(0.0, 1.0, 1); observe!(summary_b, 0.8)
  merge!(summary_a, summary_b)
  @test summary_a.counts[1] == 2
  @test GBGStreamingStatistics.Expectations._sample_mean(summary_a) ≈ 0.5
  @test GBGStreamingStatistics.Expectations._sample_variance(summary_a) ≈ ((0.2 - 0.5)^2 + (0.8 - 0.5)^2) / 2

  # A still-unimplemented accumulator type falls through to the generic erroring merge!.
  @test_throws ErrorException merge!(_DummyAccumulator(), _DummyAccumulator())
end
