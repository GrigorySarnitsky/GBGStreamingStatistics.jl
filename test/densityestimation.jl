using Test
using GBGStreamingStatistics

@testset "DensityEstimation (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:density,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
