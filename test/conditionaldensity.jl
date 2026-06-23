using Test
using GBGStreamingStatistics

@testset "ConditionalDensity (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:conditional_density, :transition_density,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
