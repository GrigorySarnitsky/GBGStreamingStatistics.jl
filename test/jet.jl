using Test
using JET
using GBGStreamingStatistics
using GBGStreamingStatistics: Histograms, MarginalDensity, Bandwidth, Regression, IncrementMoments,
                              ConditionalDensity

@testset "JET" begin
  # Static analysis of the public operations on concrete CPU summaries, with reports filtered to the
  # package's own code. Concrete summaries ⇒ `get_backend` is concretely `CPU()`, so the GPU branch of
  # the scatter kernel (only reachable with a GPU backend loaded — exercised in test/gpu/) is not part
  # of the analyzed call; the KernelAbstractions/Atomix internals are excluded by `target_modules`.
  mods = (GBGStreamingStatistics, Histograms, MarginalDensity, Bandwidth, Regression, IncrementMoments,
          ConditionalDensity)
  h = Histogram1D(-5.0, 5.0, 100; order = 4)
  h2 = Histogram1D(-5.0, 5.0, 100; order = 4)
  samples = randn(1000)

  @test_call target_modules = mods observe!(h, samples)
  @test_call target_modules = mods observe!(h, 0.3)
  @test_call target_modules = mods merge!(h, h2)
  @test_call target_modules = mods bandwidth(h)
  @test_call target_modules = mods density(h, 0.0; h = 0.5)
  # the density knobs and asymmetric kernels (closures in the kernel/start fields are the type-stability
  # risk these guard) and the bandwidth-selector ladder
  @test_call target_modules = mods density(h, 0.0; h = 0.5, kernel = gamma_kernel)
  @test_call target_modules = mods density(h, 0.0; h = 0.5, start = normal_start)
  @test_call target_modules = mods density(h, 0.0; h = 0.5, normalize = true)
  @test_call target_modules = mods density_standard_error(h, 0.0; h = 0.5)
  @test_call target_modules = mods bandwidth(h; rule = :sj)
  @test_call target_modules = mods bandwidth(h; rule = :isj)
  @test_call target_modules = mods bandwidth(h; rule = :aic)
  @test_call target_modules = mods bandwidth(h; rule = :plugin, kernel = gamma_kernel)

  g = Histogram2D(-5.0, 5.0, 50, -5.0, 5.0, 50)
  g2 = Histogram2D(-5.0, 5.0, 50, -5.0, 5.0, 50)
  xs = randn(1000); ys = randn(1000)
  @test_call target_modules = mods observe_pairs!(g, xs, ys)
  @test_call target_modules = mods merge!(g, g2)
  @test_call target_modules = mods conditional_mean(g, 0.0; h = 0.5)
  @test_call target_modules = mods conditional_variance(g, 0.0; h = 0.5)
  @test_call target_modules = mods bandwidth(g, 0.0; rule = :local)   # the response-aware local h⋆(x)

  glag = Histogram2D(-5.0, 5.0, 50, -5.0, 5.0, 50; lag = 1)
  observe_pairs!(glag, xs, ys)
  @test_call target_modules = mods increment_moment(glag, 0.0, 1; h = 0.5)
  @test_call target_modules = mods drift(glag, 0.0; h = 0.5)
  @test_call target_modules = mods diffusion(glag, 0.0; h = 0.5)
  @test_call target_modules = mods conditional_density(glag, 0.0, 0.0; hx = 0.5, hy = 0.5)
  @test_call target_modules = mods transition_density(glag, 0.0, 0.0; hx = 0.5, hy = 0.5)
end
