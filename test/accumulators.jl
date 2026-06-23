using Test
using GBGStreamingStatistics

@testset "Accumulators (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:AbstractAccumulator, :observe!, :merge!,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
