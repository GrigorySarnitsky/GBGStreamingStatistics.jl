# Backend-agnostic GPU test suite for GBGStreamingStatistics.
#
# Included by a backend entry script (test/gpu/<backend>/runtests.jl) that must define, before the
# include:
#   to_device(host_array)  -> a device array (e.g. ROCArray / CuArray)
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

# Device twins: move the backing arrays of a host summary onto the device (default constructor).
on_device(h::Histogram1D) =
  Histogram1D(h.lo, h.hi, h.scale, to_device(h.counts), to_device(h.powersums))
on_device(h::Histogram2D) =
  Histogram2D(h.xlo, h.xhi, h.xscale, h.ylo, h.yhi, h.yscale, h.nx, h.ny, h.lag,
              to_device(h.counts), to_device(h.powersums))

function run_gpu_suite(::Type{T}) where {T<:AbstractFloat}
  rtol = 1e-8                                  # power-sums are accumulated in Float64 either way
  Random.seed!(1234)

  @testset "Histogram1D observe! + merge! parity ($T)" begin
    a = randn(T, 200_000)
    b = randn(T, 200_000)
    hh = Histogram1D(-6.0, 6.0, 500; order = 5); observe!(hh, a)
    hd = on_device(Histogram1D(-6.0, 6.0, 500; order = 5)); observe!(hd, to_device(a))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)

    hh2 = Histogram1D(-6.0, 6.0, 500; order = 5); observe!(hh2, b)
    hd2 = on_device(Histogram1D(-6.0, 6.0, 500; order = 5)); observe!(hd2, to_device(b))
    merge!(hh, hh2); merge!(hd, hd2)
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "heavy-contention scatter ($T)" begin
    # Thousands of samples crammed into a few bins ⇒ the device atomics take many concurrent adds to
    # the same cell; a race or a wrong atomic would make the device differ from the CPU.
    crowded = T(0.05) .* randn(T, 50_000)
    hh = Histogram1D(-1.0, 1.0, 20; order = 4); observe!(hh, crowded)
    hd = on_device(Histogram1D(-1.0, 1.0, 20; order = 4)); observe!(hd, to_device(crowded))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
  end

  @testset "Histogram2D observe_pairs! parity ($T)" begin
    X = randn(T, 200_000); Y = randn(T, 200_000)
    hh = Histogram2D(-6.0, 6.0, 60, -6.0, 6.0, 60); observe_pairs!(hh, X, Y)
    hd = on_device(Histogram2D(-6.0, 6.0, 60, -6.0, 6.0, 60)); observe_pairs!(hd, to_device(X), to_device(Y))
    @test Array(hd.counts) == hh.counts
    @test isapprox(Array(hd.powersums), hh.powersums; rtol = rtol)
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
