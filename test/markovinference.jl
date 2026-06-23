using Test
using GBGStreamingStatistics

@testset "MarkovInference (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:increment_moment, :kramers_moyal, :drift, :diffusion, :pawula_quotient, :markov_test,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
