using Test
using Aqua
using GBGStreamingStatistics

# Package-quality checks: ambiguities, unbound type parameters, undefined exports, stale/uncompat
# dependencies, type piracy, persistent tasks. These run for real against the loadable scaffold.
@testset "Aqua" begin
  Aqua.test_all(GBGStreamingStatistics)
end
