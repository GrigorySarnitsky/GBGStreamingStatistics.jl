# Getting started

This tour shows the intended **samples in, statistics out** workflow. The architectural spine is

> empirical measure → mergeable summary → post-hoc-smoothed estimator,

and the API follows it directly.

!!! note "Status"
    The carrier and the marginal density below are implemented and runnable today; the two-time
    (conditional / transition) path at the end is forthcoming. See the
    [Architecture](architecture.md) for the design and roadmap.

## The workflow

**1. Allocate a carrier.** A [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) is a finely binned empirical measure of a
conditioning variable ``X`` — the canonical summary. It is fixed-size and lives for the whole run.

```julia
using GBGStreamingStatistics

summary = Histogram1D(-5.0, 5.0, 1000)   # fixed range and bin count; never grows with the stream
```

**2. Observe the stream.** Each batch is folded in once with [`observe!`](@ref GBGStreamingStatistics.Accumulators.observe!) and then discarded —
the samples are never stored. Accumulation is *response-agnostic*: it commits to no bandwidth, kernel,
or channel.

```julia
for batch in stream          # the consumer produces the batches; the package never holds the stream
    observe!(summary, batch)
end
```

Summaries from independent chunks or devices combine with [`merge!`](@ref Base.merge!) — the mergeable-monoid
operation that makes multi-threaded and multi-device accumulation free:

```julia
merge!(summary, other_summary)
```

**3. Estimate, post hoc.** Every reported quantity is a smoothing functional of the one stored
summary — computed after accumulation, at any bandwidth, for any channel, with no second pass:

```julia
bw = bandwidth(summary)              # a global plug-in bandwidth from the recovered σ̂
fX = density(summary, 0.0; h = bw)   # the marginal density f_X at x = 0
```

For two-time quantities (conditional and transition densities, conditional increment moments) —
*forthcoming* — the consumer forms lag-labelled pairs ``\bigl(X(t), Y(t+\tau)\bigr)`` and bins them with
[`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!) into a [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D); the snapshot buffer and lag policy stay
consumer-side.

## Where to go next

- The [Architecture](architecture.md) explains *why* this split (response-agnostic accumulation,
  post-hoc smoothing) is the load-bearing decision.
- The [API reference](api.md) documents every estimand — the population quantity before the estimator.
