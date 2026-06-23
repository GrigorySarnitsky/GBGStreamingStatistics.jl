using Test
using GBGStreamingStatistics

@testset "mergeability (monoid laws)" begin
  # Scaffold stage: every carrier is a mergeable monoid (invariant 3). Once accumulation lands these
  # assert accumulate-then-merge == accumulate-together (to round-off) and that merge! is associative.
  @test_skip accumulate_then_merge_equals_accumulate_together()
  @test_skip merge_is_associative()
end
