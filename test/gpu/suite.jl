# Backend-agnostic GPU test suite for GBGStreamingStatistics.
#
# Included by a backend entry script (test/gpu/<backend>/runtests.jl) that must define, before the
# include:
#   to_device(host_array)  -> a device array (e.g. ROCArray / CuArray)
#   DEVICE_FAMILY          -> the device array family (e.g. ROCArray), the summaries' first constructor argument
#   device_functional()    -> Bool, whether a usable device is present
#   GPU_BACKEND::String     -> a human-readable backend label
# and must already have called <Backend>.allowscalar(false). With scalar indexing disabled, any
# scalar index into device memory inside a *package method* throws, so the suite mechanically enforces
# the no-scalar-indexing invariant (CLAUDE.md §4) rather than trusting it by inspection.
#
# The same `observe!` / `merge!` / `observe_pairs!` are the *one portable kernel* run on the CPU
# (the trusted reference) and on the device; the suite asserts the device result equals the CPU one
# (counts exactly — integer atomics are order-independent; power-sums to round-off — float atomic-adds
# accumulate in a nondeterministic order).

using Test
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms
using Random
const H = Histograms


function run_gpu_suite(::Type{T}) where {T<:AbstractFloat}
  rtol = 1e-8                                  # power-sums are accumulated in Float64 either way
  Random.seed!(1234)

  @testset "Histogram1D observe! + merge! parity ($T)" begin
    a = randn(T, 200_000)
    b = randn(T, 200_000)
    hh = Histogram1D(-6.0, 6.0, 500; order = 5); observe!(hh, a)
    hd = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 500; order = 5); observe!(hd, to_device(a))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)

    hh2 = Histogram1D(-6.0, 6.0, 500; order = 5); observe!(hh2, b)
    hd2 = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 500; order = 5); observe!(hd2, to_device(b))
    merge!(hh, hh2); merge!(hd, hd2)
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "cone parity — Float-count observe!, scale!, weighted merge! ($T)" begin
    # the Float-count variant routes observe! through the GPU atomic FLOAT add (distinct from the
    # integer-count atomic path); scale! and the weighted merge! are device broadcasts.
    a = randn(T, 200_000); b = randn(T, 200_000)
    hh = Histogram1D(-6.0, 6.0, 500; order = 5, count_type = Float64); observe!(hh, a)
    hd = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 500; order = 5, count_type = Float64); observe!(hd, to_device(a))
    @test isapprox(Array(hd.counts), hh.counts; rtol = rtol)
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)

    scale!(hh, 0.7); scale!(hd, 0.7)
    @test isapprox(Array(hd.counts), hh.counts; rtol = rtol)
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)

    gh = Histogram1D(-6.0, 6.0, 500; order = 5, count_type = Float64); observe!(gh, b)
    gd = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 500; order = 5, count_type = Float64); observe!(gd, to_device(b))
    merge!(hh, gh; weight = 0.3); merge!(hd, gd; weight = 0.3)
    @test isapprox(Array(hd.counts), hh.counts; rtol = rtol)
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "heavy-contention scatter ($T)" begin
    # Thousands of samples crammed into a few bins ⇒ the device atomics take many concurrent adds to
    # the same cell; a race or a wrong atomic would make the device differ from the CPU.
    crowded = T(0.05) .* randn(T, 50_000)
    hh = Histogram1D(-1.0, 1.0, 20; order = 4); observe!(hh, crowded)
    hd = Histogram1D(DEVICE_FAMILY, -1.0, 1.0, 20; order = 4); observe!(hd, to_device(crowded))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "Histogram2D observe_pairs! parity ($T)" begin
    X = randn(T, 200_000); Y = randn(T, 200_000)
    hh = Histogram2D(-6.0, 6.0, 60, -6.0, 6.0, 60); observe_pairs!(hh, X, Y)
    hd = Histogram2D(DEVICE_FAMILY, -6.0, 6.0, 60, -6.0, 6.0, 60); observe_pairs!(hd, to_device(X), to_device(Y))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "a field-shaped device batch folds in place ($T)" begin
    # A consumer's field is a 3-D device array (or a view of one); observe! reads it where it lies.
    field = randn(T, 32, 32, 16)
    hh = Histogram1D(-6.0, 6.0, 400); observe!(hh, field)
    hd = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 400); observe!(hd, to_device(field))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
    hv = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 400); observe!(hv, view(to_device(field), :, :, 1))
    hhv = Histogram1D(-6.0, 6.0, 400); observe!(hhv, field[:, :, 1])
    @test Array(hv.counts) == hhv.counts
  end

  @testset "readouts of a device summary equal the host summary's ($T)" begin
    x = randn(T, 100_000); y = T(0.6) .* x .+ T(0.8) .* randn(T, 100_000)
    hh = Histogram1D(-6.0, 6.0, 600); observe!(hh, x)
    hd = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 600); observe!(hd, to_device(x))
    xs = collect(range(-2.0, 2.0; length = 9))
    @test bandwidth(hd) ≈ bandwidth(hh)
    @test density(hd, xs; h = 0.3) ≈ density(hh, xs; h = 0.3)
    @test density(hd, xs; h = 0.3, degree = 1) ≈ density(hh, xs; h = 0.3, degree = 1)
    @test density_standard_error(hd, xs; h = 0.3) ≈ density_standard_error(hh, xs; h = 0.3)
    @test mean(hd) ≈ mean(hh)
    @test var(hd) ≈ var(hh)
    @test expectation(hd, cos) ≈ expectation(hh, cos)
    jh = Histogram2D(-6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(jh, x, y)
    jd = Histogram2D(DEVICE_FAMILY, -6.0, 6.0, 80, -6.0, 6.0, 80); observe_pairs!(jd, to_device(x), to_device(y))
    @test conditional_mean(jd, xs; h = 0.4) ≈ conditional_mean(jh, xs; h = 0.4)
    @test conditional_variance(jd, xs; h = 0.4) ≈ conditional_variance(jh, xs; h = 0.4)
    @test conditional_density(jd, 0.3, xs; h_x = 0.4, h_y = 0.3) ≈ conditional_density(jh, 0.3, xs; h_x = 0.4, h_y = 0.3)
  end
end

@testset "allocating constructors build on the device — $GPU_BACKEND" begin
  if device_functional()
    h = Histogram1D(DEVICE_FAMILY, -1.0, 1.0, 50; order = 6, count_type = Float64, accumulator_type = Float32)
    @test h.counts isa DEVICE_FAMILY && h.powersums isa DEVICE_FAMILY
    @test eltype(h.counts) === Float64 && eltype(h.powersums) === Float32
    @test size(h.powersums) == (6, 50)
    @test all(iszero, Array(h.counts)) && all(iszero, Array(h.powersums))
    host = Histogram1D(-1.0, 1.0, 50; order = 6)
    d = Histogram1D(DEVICE_FAMILY, -1.0, 1.0, 50; order = 6)
    @test (d.lower, d.upper, d.scale) == (host.lower, host.upper, host.scale)
    j = Histogram2D(DEVICE_FAMILY, -1.0, 1.0, 7, -2.0, 2.0, 9; lag = 3)
    @test j.counts isa DEVICE_FAMILY && length(j.counts) == 63 && size(j.powersums) == (5, 63)
    @test j.lag == 3 && all(iszero, Array(j.counts))
    @test_throws ArgumentError Histogram1D(DEVICE_FAMILY{Float32}, -1.0, 1.0, 50)
    @test_throws ArgumentError Histogram2D(DEVICE_FAMILY{Float64}, -1.0, 1.0, 7, -1.0, 1.0, 7)
  end
end

@testset "a summary and its data on different backends are refused — $GPU_BACKEND" begin
  if device_functional()
    a = randn(1000)
    host = Histogram1D(-6.0, 6.0, 100); dev = Histogram1D(DEVICE_FAMILY, -6.0, 6.0, 100)
    @test_throws ArgumentError observe!(host, to_device(a))
    @test_throws ArgumentError observe!(dev, a)
    @test_throws ArgumentError merge!(dev, host)
    @test_throws ArgumentError merge!(host, dev)
    @test all(iszero, host.counts) && all(iszero, Array(dev.counts))   # a refusal folds nothing
    jh = Histogram2D(-6.0, 6.0, 20, -6.0, 6.0, 20); jd = Histogram2D(DEVICE_FAMILY, -6.0, 6.0, 20, -6.0, 6.0, 20)
    @test_throws ArgumentError observe_pairs!(jh, to_device(a), to_device(a))
    @test_throws ArgumentError observe_pairs!(jd, a, a)
    @test_throws ArgumentError observe_pairs!(jd, to_device(a), a)
    @test_throws ArgumentError merge!(jd, jh)
  end
end

@testset "GPU parity — $GPU_BACKEND" begin
  if !device_functional()
    @warn "no functional $GPU_BACKEND device; skipping the parity suite"
  else
    for T in (Float64, Float32)
      run_gpu_suite(T)
    end
  end
end
