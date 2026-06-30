# Accumulators — the foundation (universal lower half): the mergeable-monoid accumulator interface.
# The concrete summaries and their KernelAbstractions + Atomix scatter kernels live in `Histograms`.
# Conceptual model and conventions: `docs/ARCHITECTURE.md` (rationale), `CLAUDE.md` (operational rules).

module Accumulators

import Base: merge!

export AbstractAccumulator, observe!, merge!, scale!


"""
    AbstractAccumulator

Supertype of every **mergeable summary** ``T(\\hat\\mu_n)`` of the [empirical measure](https://en.wikipedia.org/wiki/Empirical_measure)
``\\hat\\mu_n = \\frac1n \\sum_{i=1}^{n} \\delta_{X_i}`` of an observed stream ``X_1, X_2, \\dots``
(``\\delta_x`` the unit point mass at ``x``, ``n`` the number of samples seen).

A summary is **finite-dimensional** — the ``n`` point masses are never stored — **mergeable** (an
additive monoid over the stream; see [`merge!`](@ref Base.merge!)), and **rich enough** that every
reported estimator is recoverable from it post hoc, with no second pass over the data. The summary is
accumulated online and is *estimand-agnostic*; all smoothing (kernel, bandwidth) is applied
post hoc. Concrete summaries ([`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D),
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D)) subtype this.
"""
abstract type AbstractAccumulator end


"""
    observe!(summary::AbstractAccumulator, batch) -> summary

Fold a `batch` of samples into `summary` in place, advancing the accumulated summary
``T(\\hat\\mu_n)`` of the empirical measure (see
[`AbstractAccumulator`](@ref GBGStreamingStatistics.Accumulators.AbstractAccumulator)). Accumulation
is **estimand-agnostic**: it knows nothing of which estimator will later be asked for, nor of the
bandwidth or kernel — the architectural invariant that delivers *no reruns*.

The binning is an **atomic scatter-add**: one [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl) + [Atomix](https://github.com/JuliaConcurrent/Atomix.jl) kernel — one work-item per
sample — accumulates, per bin, the **count** and the **power-sums**
``S_p = \\sum_{i\\in\\text{bin}}(X_i - c_j)^p`` about the bin center ``c_j``, directly into the
preallocated summary (no scratch). It runs identically on a CPU `Array` and any GPU backend, and is
allocation-free up to the kernel launch.
Because ``X_i - c_j`` is small within a bin, the local power-sums are well-conditioned (no compensated
summation needed); global central moments are recovered post hoc by the shift formula.
Accumulated in the summary's accumulator type (`Float64` by default),
never narrowed to the input element type — *accumulate wider than the data*.

References:

- D. W. Scott (2015), [*Multivariate Density Estimation*](https://doi.org/10.1002/9781118575574), 2nd ed.
- M. D. Cattaneo et al. (2024), [*On binscatter*](https://doi.org/10.1257/aer.20221576) (the local-moment / WARPing construction).
"""
function observe! end


"""
    merge!(into::AbstractAccumulator, from::AbstractAccumulator; weight = 1) -> into

Combine two summaries into `into`, forming the summary of the concatenated stream — the monoid
operation that makes accumulation additive. For disjoint sample blocks ``A`` and ``B``,
``T(A \\uplus B) = T(A) \\oplus T(B)``: accumulation ``T`` is a **monoid homomorphism** from sample
blocks (under disjoint union ``\\uplus``) to summaries (under `merge!`), on the **unnormalized** counts
and power-sums — the normalized ``\\hat\\mu`` would reweight, not add. Because compatible summaries
share the same binning (hence the same bin centers),
their **counts and power-sums simply add elementwise**, so `merge!` is exact and trivially parallel —
this is what makes multi-device reduction and time-chunk stitching free (combine per-device or
per-time-chunk summaries by addition).

**Weighted form — the cone.** `merge!(into, from; weight)` forms
``\\texttt{into} \\mathrel{+}= \\texttt{weight} \\cdot \\texttt{from}``; the default ``\\texttt{weight} = 1``
is the plain monoid addition above. Because the unnormalized summary is also **scalable**
([`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!)), the summaries form a **convex cone**
(closed under nonnegative linear combinations), and the weighted merge realizes the **hyperstatistical**
parameter-mixture density (superstatistics; Beck & Cohen 2003)

```math
P(x) = \\int f(\\beta)\\,P(x \\mid \\beta)\\,\\mathrm{d}\\beta \\;\\approx\\; \\sum_k w_k\\,f_k(x),
\\qquad w_k = f(\\beta_k)\\,\\Delta\\beta,
```

a superposition over a slowly varying parameter ``\\beta`` with mixing law ``f(\\beta)``: keep one summary
``T_k`` per value ``\\beta_k`` (accumulated independently), then **after** the runs choose ``f(\\beta)`` and
assemble ``\\sum_k (w_k / N_k)\\,T_k`` by weighted `merge!` (``N_k = \\texttt{sum(T\\_k.counts)}`` normalizes
each component to a probability measure; the conditional / transition estimators then carry the local
mixing weight automatically, being ratios of linear functionals). This is **post hoc** — re-weight without
re-accumulating, exactly as bandwidth is a post-hoc choice. A weighted merge (``\\texttt{weight} \\ne 1``)
requires a Float-count `into` (`count = Float64`).

`merge!` is intrinsically a finite-sample operation, so it has no population estimand of its own; its
specification is the monoid (and, weighted, the cone) law, tested as *accumulate-then-merge `=`
accumulate-together*. Extends `Base.merge!`.
"""
merge!(into::AbstractAccumulator, from::AbstractAccumulator; weight = 1) =
  error("merge! is not implemented for $(typeof(into))")


"""
    scale!(summary::AbstractAccumulator, weight) -> summary

Reweight `summary` in place by a nonnegative `weight` — multiply every count and power-sum by `weight`,
giving the valid summary of the reweighted measure ``\\texttt{weight} \\cdot \\hat\\mu``. With
[`merge!`](@ref Base.merge!) this makes the unnormalized summaries a **convex cone** (closed under
nonnegative linear combinations), not merely an additive monoid: `scale!` is the scalar action, `merge!`
the addition.

Two uses. **Exponential forgetting** of a non-stationary stream — `scale!(summary, λ)` with
``0 < \\lambda < 1`` before folding in the next block down-weights the past geometrically (an effective
window ``1/(1-\\lambda)``). **Hyperstatistics** — the scalar building block of the parameter-mixture
density assembled by weighted [`merge!`](@ref Base.merge!) (the estimand and recipe are documented there;
Beck & Cohen 2003).

Requires a **Float-count summary** (`count = Float64`; the default `UInt64` counts cannot hold a
fractional ``\\texttt{weight} \\cdot n``) — an `ArgumentError` otherwise. ``\\texttt{weight} = 1`` is the
identity. Like `merge!`, an intrinsically finite-sample operation with no population estimand; specified
by the cone laws (`scale!(scale!(s, a), b)` equals scaling by ``ab``; `scale!` distributes over `merge!`)
and tested as such.
"""
function scale! end

end # module Accumulators
