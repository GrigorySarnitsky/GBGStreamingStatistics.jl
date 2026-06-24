using Test
using JET
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms, DensityEstimation, Bandwidth

@testset "JET" begin
  # Static analysis of the public operations on concrete CPU carriers, with reports filtered to the
  # package's own code. Concrete carriers ⇒ `get_backend` is concretely `CPU()`, so the GPU branch of
  # the scatter kernel (only reachable with a GPU backend loaded — exercised in test/gpu/) is not part
  # of the analyzed call; the KernelAbstractions/Atomix internals are excluded by `target_modules`.
  mods = (GBGStreamingStatistics, Histograms, DensityEstimation, Bandwidth)
  h = Histogram1D(-5.0, 5.0, 100; order = 4)
  h2 = Histogram1D(-5.0, 5.0, 100; order = 4)
  samples = randn(1000)

  @test_call target_modules = mods observe!(h, samples)
  @test_call target_modules = mods observe!(h, 0.3)
  @test_call target_modules = mods merge!(h, h2)
  @test_call target_modules = mods bandwidth(h)
  @test_call target_modules = mods density(h, 0.0; h = 0.5)

  g = Histogram2D(-5.0, 5.0, 50, -5.0, 5.0, 50)
  g2 = Histogram2D(-5.0, 5.0, 50, -5.0, 5.0, 50)
  xs = randn(1000); ys = randn(1000)
  @test_call target_modules = mods observe_pairs!(g, xs, ys)
  @test_call target_modules = mods merge!(g, g2)
end
