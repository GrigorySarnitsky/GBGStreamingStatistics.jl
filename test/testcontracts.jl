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

export check_point_estimator, check_standard_error, check_se_bandwidth_scaling

"""
    check_point_estimator(f, merged, full; points, outside = nothing, sentinel = NaN, rtol = 1e-8)

Assert the universal query-point-estimator contract. `f(summary, x)` is the estimator with its bandwidths
fixed (a closure); `merged` is a summary `merge!`d from disjoint parts of a stream (**two or more** — pass
a 3-way merge to exercise associativity) and `full` the summary accumulated from the whole stream (so the
two hold the same data); `points` are in-support query points.

Checks: (1) **scalar == vector** — `f(full, points)` equals `[f(full, x) for x in points]`; (2)
**functional of the merged summary** — `f(merged, x) ≈ f(full, x)` (to `rtol`; this holds for both a point
estimate and its standard error, since `merged` and `full` hold the same data). If `outside` is given, also
(3) the **off-support sentinel** — `f(full, outside)` is `NaN` (default) or `≈ sentinel`.
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


"""
    check_standard_error(se, s; points, K = 4)

Assert the **standard-error contract**. `se(summary, x)` is a standard error with its options fixed (a
closure); `s` an in-support summary; `points` in-support query points. Checks: **positivity** (`se > 0` at
each point); **scalar == vector**; and the **`n^{-1/2}` scaling** — merging `K` copies of `s` (`K×` the
samples, the same shape) shrinks the standard error by exactly `√K` (`se(K·s) ≈ se(s)/√K`), the
deterministic form of `SE ∝ n^{-1/2}`.
"""
function check_standard_error(se, s; points, K = 4)
  @testset "standard-error contract" begin
    @test all(x -> se(s, x) > 0, points)
    @test se(s, points) ≈ [se(s, x) for x in points]
    sK = deepcopy(s)
    for _ in 2:K
      merge!(sK, s)
    end
    for x in points
      @test se(sK, x) ≈ se(s, x) / sqrt(K)  rtol = 1e-6
    end
  end
end

"""
    check_se_bandwidth_scaling(se_h, s, x; h, rtol = 0.06)

Assert the **``h^{-1/2}`` bandwidth scaling** of a *kernel* standard error. `se_h(summary, x, h)` is the SE
at bandwidth `h`. Unlike the exact ``n^{-1/2}`` scaling, this is **asymptotic** — changing `h` also moves
`f̂` and the second variance term `−f̂²` (weighted `1/n`, or locally for a conditional estimate) carries no
`1/h` — so it holds *at small `h` on a smooth region*, where that term is negligible (pass a flat/uniform `s`
and an interior `x`, or a point of low density): doubling the bandwidth shrinks the leading-variance SE by
``\\sqrt2`` (`se_h(s, x, 2h) ≈ se_h(s, x, h)/√2`).
"""
function check_se_bandwidth_scaling(se_h, s, x; h, rtol = 0.06)
  @testset "SE bandwidth (h^{-1/2}) scaling" begin
    @test se_h(s, x, 2 * h) ≈ se_h(s, x, h) / sqrt(2)  rtol = rtol
  end
end

end # module TestContracts
