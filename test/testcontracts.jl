# The universal contract every query-point estimator must satisfy — applied by *construction* (one call),
# not by convention, so a new estimator cannot silently omit a check (the failure mode that left
# `conditional_mean` without a mergeability test in the first place). A query-point estimator is a
# functional of the summary evaluated at a query point; the universal contract is: its vector form equals
# its scalar form, it is a deterministic functional of the *merged* summary, and — where the estimand is
# undefined off the data's support — it returns the agreed sentinel. The *estimand-specific*
# characteristics (an independent analytic oracle, normalization to 1 for a density, the Kramers–Moyal
# identities, the convergence rate) are layered on top, per estimator (see docs/TESTING.md).
module TestContracts

using Test

export check_point_estimator

"""
    check_point_estimator(f, merged, full; points, outside = nothing, sentinel = NaN, rtol = 1e-8)

Assert the universal query-point-estimator contract. `f(summary, x)` is the estimator with its bandwidths
fixed (a closure); `merged` is a summary `merge!`d from disjoint halves of a stream and `full` the summary
accumulated from the whole stream (so the two hold the same data); `points` are in-support query points.

Checks: (1) **scalar == vector** — `f(full, points)` equals `[f(full, x) for x in points]`; (2)
**functional of the merged summary** — `f(merged, x) ≈ f(full, x)` (to `rtol`). If `outside` is given,
also (3) the **off-support sentinel** — `f(full, outside)` is `NaN` (default) or `≈ sentinel`.
"""
function check_point_estimator(f, merged, full; points, outside = nothing, sentinel = NaN, rtol = 1e-8)
  @testset "universal estimator contract" begin
    @test f(full, points) ≈ [f(full, x) for x in points]              # scalar / vector consistency
    for x in points
      @test f(merged, x) ≈ f(full, x)  rtol = rtol                    # deterministic functional of merge
    end
    if outside !== nothing
      if sentinel isa Real && isnan(sentinel)
        @test isnan(f(full, outside))
      else
        @test f(full, outside) ≈ sentinel
      end
    end
  end
end

end # module TestContracts
