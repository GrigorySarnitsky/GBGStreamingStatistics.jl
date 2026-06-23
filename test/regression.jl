using Test
using GBGStreamingStatistics

@testset "Regression (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:conditional_mean, :conditional_variance,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
