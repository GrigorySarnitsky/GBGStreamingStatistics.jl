using Test
using GBGStreamingStatistics

@testset "Histograms (scaffold)" begin
  # Scaffold stage: the public surface is declared; estimand-first oracle, rate, and mergeability
  # assertions land with the implementation (docs/TESTING.md). For now, assert the surface exists.
  for name in (:Histogram1D, :Histogram2D, :bin_index, :observe_pairs!,)
    @test isdefined(GBGStreamingStatistics, name)
  end
end
