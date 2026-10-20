using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random, Statistics
const H = Histograms

# A stand-in for a device array on a machine without one: a vector whose backend is not the CPU's, so the
# backend check can be exercised here (the real device cases run in the gated GPU suite, test/gpu/).
struct ElsewhereVector <: AbstractVector{Float64}
  data::Vector{Float64}
end
Base.size(v::ElsewhereVector) = size(v.data)
Base.getindex(v::ElsewhereVector, i::Int) = v.data[i]
struct ElsewhereBackend end
H.get_backend(::ElsewhereVector) = ElsewhereBackend()

@testset "Histogram1D" begin
  @testset "construction & binning" begin
    h = Histogram1D(-5.0, 5.0, 10)
    @test H._nbins(h) == 10
    @test H._order(h) == 4
    @test_throws ArgumentError Histogram1D(5.0, -5.0, 10)
    @test_throws ArgumentError Histogram1D(-5.0, 5.0, 0)
    @test_throws ArgumentError Histogram1D(-5.0, 5.0, 10; order = 1)
  end

  @testset "the array-family constructor" begin
    # `Histogram1D(Array, …)` is the host form spelled out: the same type, fields and zeroed storage.
    for kw in ((;), (; order = 6), (; count_type = Float64, accumulator_type = Float32))
      a = Histogram1D(-2.0, 3.0, 7; kw...)
      b = Histogram1D(Array, -2.0, 3.0, 7; kw...)
      @test typeof(a) == typeof(b)
      @test (a.lower, a.upper, a.scale) == (b.lower, b.upper, b.scale)
      @test a.counts == b.counts && a.powersums == b.powersums
      @test all(iszero, b.counts) && all(iszero, b.powersums)
    end
    @test_throws ArgumentError Histogram1D(Array{Float64}, -1.0, 1.0, 4)   # element types are keywords
    @test_throws ArgumentError Histogram1D(Vector{UInt64}, -1.0, 1.0, 4)
    @test_throws ArgumentError Histogram1D(AbstractArray, -1.0, 1.0, 4)    # not an allocatable family
    @test_throws ArgumentError Histogram1D(Array, 1.0, -1.0, 4)            # the range checks still apply
    @test_throws ArgumentError Histogram1D(Array, -1.0, 1.0, 0)
  end

  @testset "a batch on another backend is refused" begin
    h = Histogram1D(-5.0, 5.0, 10)
    @test_throws ArgumentError observe!(h, ElsewhereVector([0.5, 1.5]))
    @test all(iszero, h.counts)                       # a refusal folds nothing
    # host batches without a KernelAbstractions backend (a range, a view) are host data and accepted
    observe!(h, 1:3); observe!(h, view([0.1, 0.2, 0.3], 1:2)); observe!(h, reshape([0.1, 0.2], 1, 2))
    @test H._total(h) == 7
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

    # Exact agreement with the direct sample central moments: the power-sums about the bin center recover
    # the moments by the shift formula, which is exact arithmetic (no approximation in the reduction).
    m = mean(data)
    direct(p) = sum((data .- m) .^ p) / length(data)
    for p in 2:6
      @test GBGStreamingStatistics.Expectations._sample_central_moment(h, p) ≈ direct(p) rtol = 1e-8 atol = 1e-12
    end

    # Agreement with the true Gaussian moments, to sampling error.
    @test GBGStreamingStatistics.Expectations._sample_mean(h) ≈ 0 atol = 5e-3
    @test GBGStreamingStatistics.Expectations._sample_variance(h) ≈ 1 rtol = 1e-2
    @test abs(GBGStreamingStatistics.Expectations._sample_central_moment(h, 3)) < 1e-2          # skewness ≈ 0
    @test GBGStreamingStatistics.Expectations._sample_central_moment(h, 4) ≈ 3 rtol = 2e-2      # kurtosis
    @test abs(GBGStreamingStatistics.Expectations._sample_central_moment(h, 5)) < 5e-2          # ≈ 0
    @test GBGStreamingStatistics.Expectations._sample_central_moment(h, 6) ≈ 15 rtol = 5e-2
  end

  @testset "Float32 input accumulates in Float64" begin
    Random.seed!(3)
    h = Histogram1D(-8.0, 8.0, 1000; order = 4)
    observe!(h, randn(Float32, 300_000))
    @test eltype(h.powersums) === Float64
    @test eltype(h.counts) === UInt64
    @test GBGStreamingStatistics.Expectations._sample_variance(h) ≈ 1 rtol = 3e-2
  end

  @testset "public moment accessors (mean/var/std/central_moment)" begin
    Random.seed!(11)
    data = randn(1_000_000)
    h = Histogram1D(-8.0, 8.0, 1600; order = 6)
    observe!(h, data)
    @test H._total(h) == length(data)                           # all samples in range
    @test mean(h) == GBGStreamingStatistics.Expectations._sample_mean(h)                          # public == internal
    @test central_moment(h, 4) == GBGStreamingStatistics.Expectations._sample_central_moment(h, 4)
    # var/std default to Bessel's (n-1) correction, matching Statistics on the same data
    @test var(h) ≈ var(data) rtol = 1e-6                        # corrected == Statistics.var
    @test std(h) ≈ std(data) rtol = 1e-6                        # corrected == Statistics.std
    @test var(h; corrected = false) == GBGStreamingStatistics.Expectations._sample_variance(h)    # 1/n == internal helper ...
    @test var(h; corrected = false) == central_moment(h, 2)     # ... == the order-2 central moment
    @test std(h; corrected = false) == GBGStreamingStatistics.Expectations._sample_std(h)
    @test var(h) > var(h; corrected = false)                    # n/(n-1) > 1
    @test mean(h) ≈ 0 atol = 5e-3                               # Gaussian oracle
    @test var(h) ≈ 1 rtol = 2e-2
    @test std(h) ≈ 1 rtol = 2e-2
    @test central_moment(h, 4) ≈ 3 rtol = 3e-2
    @test_throws ArgumentError central_moment(h, 7)             # p > order
  end
end
