using Test
using GBGStreamingStatistics

@testset "allocation audit (@allocated)" begin
  # Scaffold stage: the hot path is allocation-free (invariant 4). Once accumulation lands this
  # asserts, via @allocated on the CPU, that observe!/observe_pairs! allocate 0 bytes on a
  # preallocated summary — the contract that keeps a long GPU run from wedging.
  @test_skip observe_is_allocation_free()
  @test_skip observe_pairs_is_allocation_free()
end
