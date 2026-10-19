# Accumulators — the foundation (universal lower half): the mergeable-monoid accumulator interface.
# The concrete summaries and their KernelAbstractions + Atomix scatter kernels live in `Histograms`.
# Conceptual model and conventions: `docs/ARCHITECTURE.md` (rationale), `CLAUDE.md` (operational rules).

module Accumulators

import Base: merge!

export AbstractAccumulator, observe!, merge!, scale!


"""
    AbstractAccumulator

It is the supertype of every mergeable summary ``T(n\\mathbb{P}_n)`` of the [empirical measure](https://en.wikipedia.org/wiki/Empirical_measure)
``\\mathbb{P}_n`` of an observed stream ``X_1, X_2, \\dots``, which gives a set ``A`` the share of samples in it,
``\\mathbb{P}_n(A) = \\#\\{i : X_i \\in A\\}/n`` (``n`` the number of samples seen).

A summary is finite-dimensional — the ``n`` samples are never stored — mergeable and
scalable (a convex cone over the stream; see [`merge!`](@ref Base.merge!) and [`scale!`](@ref)), and
rich enough that every reported estimator is recoverable from it post hoc, with no second pass over the
data. The summary is accumulated online and is *estimand-agnostic*, independent of what is later
estimated; all smoothing (kernel, bandwidth) is applied post hoc. Concrete summaries ([`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D),
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D)) subtype this.
"""
abstract type AbstractAccumulator end


"""
    observe!(summary::AbstractAccumulator, batch) -> summary

Fold a `batch` of samples (a vector, or on a CPU summary a single sample) into `summary` in place, advancing the accumulated summary
``T(n\\mathbb{P}_n)`` of the empirical measure (see
[`AbstractAccumulator`](@ref GBGStreamingStatistics.Accumulators.AbstractAccumulator)). Accumulation
is estimand-agnostic: it knows nothing of which estimator will later be asked for, nor of the
bandwidth or smoothing kernel, so nothing is recomputed from the raw stream.

The binning is an atomic scatter-add: one [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl) + [Atomix](https://github.com/JuliaConcurrent/Atomix.jl) compute kernel — one work-item per
sample — accumulates per bin the count and the power-sums ``S_p`` of the samples' offsets from the bin center
that a [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) carries, directly into the
preallocated summary. The batch form runs identically on a CPU `Array` and any GPU backend, is
allocation-free up to the kernel launch, and accumulates *wider than the data* (in the element type of the
power-sums, `Float64` by default, never narrowed to the input element type).

# References

- D. W. Scott (2015), [*Multivariate Density Estimation*](https://doi.org/10.1002/9781118575574), 2nd ed.
- M. D. Cattaneo et al. (2024), [*On binscatter*](https://doi.org/10.1257/aer.20221576) (the local-moment construction, WARPing — weighted averaging of rounded points).
"""
function observe! end


"""
    merge!(into::AbstractAccumulator, from::AbstractAccumulator; weight = 1) -> into

Combine two summaries into `into`, forming the summary of the concatenated stream — the monoid
operation that makes accumulation additive. For disjoint sample blocks ``\\mathcal{A}`` and ``\\mathcal{B}``,
``T(\\mathcal{A} \\uplus \\mathcal{B}) = T(\\mathcal{A}) \\oplus T(\\mathcal{B})``: accumulation ``T``, which maps a block to
its summary, is a *monoid homomorphism* from sample blocks (under disjoint union ``\\uplus``) to summaries
(under `merge!`, the circled plus), on the unnormalized counts and power-sums — the normalized empirical
measure ``\\mathbb{P}_n`` would reweight, not add. Because compatible summaries share the same binning
(hence the same bin centers), their counts and power-sums simply add elementwise, so `merge!` is exact
and trivially parallel — this is what makes multi-device reduction and time-chunk stitching free.

**The weighted form is the cone.** `merge!(into, from; weight)`, the multiplier `weight` written ``\\lambda``, forms
``\\texttt{into} \\mathrel{+}= \\lambda \\cdot \\texttt{from}``; the default ``\\lambda = 1``
is the plain monoid addition above. Because the unnormalized summary is also scalable
([`scale!`](@ref GBGStreamingStatistics.Accumulators.scale!)), the summaries form a *convex cone*
(closed under nonnegative linear combinations), and the weighted merge realizes the hyperstatistical
parameter-mixture density (superstatistics; Beck & Cohen 2003)

```math
f_X(x) = \\int \\pi(\\theta) f_\\theta(x)\\,\\mathrm{d}\\theta \\approx \\sum_k \\pi_k f_{\\theta_k}(x),
\\qquad \\pi_k = \\pi(\\theta_k) (\\theta_{k+1} - \\theta_k),
```

a superposition, over a slowly varying parameter ``\\theta`` with mixing law ``\\pi(\\theta)``, of the densities of
the samples drawn at each fixed parameter value: keep one summary
``T_k`` per value ``\\theta_k``, then, after the runs, choose ``\\pi(\\theta)`` and
assemble ``\\sum_k (\\pi_k / n(T_k)) T_k`` by weighted `merge!` (the total count ``n(T_k) = \\texttt{sum(T\\_k.counts)}`` normalizes
each component to a probability measure; the conditional and transition estimators, read off a
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D), then weight each component by its
mass near the conditioning value, being ratios of linear functionals). This is post hoc, like the
bandwidth. A weighted merge (``\\lambda \\ne 1``) requires `into` to have floating-point counts
(built with `count_type = Float64`).

Being intrinsically a finite-sample operation, `merge!` has no population estimand of its own; its
specification is the monoid (and, weighted, the cone) law, tested as *accumulate-then-merge equals
accumulate-together*. Extends `Base.merge!`.

# Arguments

- `into` — the summary that receives the merge, updated in place.
- `from` — the summary folded into `into`; both share one binning (range, bin count, power-sum order), and `from`
  is left unchanged.
- `weight` — the nonnegative weight of `from` in the cone combination (default `1`, the plain monoid
  addition); a weight other than `1` needs floating-point counts in `into`.
"""
merge!(into::AbstractAccumulator, from::AbstractAccumulator; weight = 1) =
  error("merge! is not implemented for $(typeof(into))")


"""
    scale!(summary::AbstractAccumulator, weight) -> summary

Reweight `summary` in place by a nonnegative multiplier `weight`, written ``\\lambda`` — multiply every count and power-sum by it,
giving the summary of the reweighted counting measure ``\\lambda \\cdot n\\mathbb{P}_n`` of the samples
seen. With
[`merge!`](@ref Base.merge!) this makes the unnormalized summaries a convex cone (closed under
nonnegative linear combinations), not merely an additive monoid: `scale!` is the scalar action, `merge!`
the addition.

It has two uses:

- **Exponential forgetting** of a non-stationary stream: `scale!(summary, λ)` with
  ``0 < \\lambda < 1`` before folding in the next block down-weights the past geometrically (an effective
  window ``1/(1-\\lambda)``).
- **Hyperstatistics** uses it as the scalar building block of the
  parameter-mixture density assembled by weighted `merge!` (the estimand and recipe are documented there;
  Beck & Cohen 2003).

It requires floating-point counts (`count_type = Float64`; the default `UInt64` counts cannot hold a
fractional count ``\\lambda n``) and throws an `ArgumentError` otherwise. The value ``\\lambda = 1``
is the identity. Like `merge!`, it is an intrinsically finite-sample operation with no population estimand,
specified by the cone laws (scaling by ``\\lambda`` then by ``\\lambda'`` equals scaling by
``\\lambda\\lambda'``; `scale!` distributes over `merge!`).
"""
function scale! end

end # module Accumulators
