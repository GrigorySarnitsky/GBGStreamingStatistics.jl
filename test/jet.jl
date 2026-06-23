using Test
using JET
using GBGStreamingStatistics

@testset "JET" begin
  # Static analysis of the package. The scaffold has few method bodies, so this mainly guards the
  # loadable surface; per-operation `@test_call` / `@test_opt` checks land with the implementations.
  JET.test_package(GBGStreamingStatistics; target_modules = (
    GBGStreamingStatistics,
    GBGStreamingStatistics.Accumulators,
    GBGStreamingStatistics.Histograms,
    GBGStreamingStatistics.DensityEstimation,
    GBGStreamingStatistics.Regression,
    GBGStreamingStatistics.ConditionalDensity,
    GBGStreamingStatistics.MarkovInference,
    GBGStreamingStatistics.Bandwidth,
  ))
end
