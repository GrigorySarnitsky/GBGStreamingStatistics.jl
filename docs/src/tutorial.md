# Getting started

This tour shows the intended **samples in, statistics out** workflow. The architectural spine is

> empirical measure → mergeable summary → post-hoc-smoothed estimator,

and the API follows it directly.

!!! note "Status"
    Everything in this tour is implemented and runnable today, single-node on CPU or GPU. The
    Kramers–Moyal / SDE-reading functionals (drift and diffusion), the higher Kramers–Moyal
    coefficients, and the Pawula / Markov diagnostics are not currently provided but are planned; see
    the [Architecture](architecture.md) for the design and roadmap.

## The workflow

**1. Allocate a summary.** A [`Histogram1D`](@ref) is a finely binned empirical measure of a
conditioning variable ``X`` — the canonical summary. It is fixed-size and lives for the whole run.

```julia
using GBGStreamingStatistics

summary = Histogram1D(-5.0, 5.0, 1000)   # fixed range and bin count; never grows with the stream
```

**2. Observe the stream.** Each batch is folded in once with [`observe!`](@ref) and then discarded —
the samples are never stored. Accumulation is *estimand-agnostic*: it commits to no channel, bandwidth,
or kernel.

```julia
for batch in stream          # the consumer produces the batches; the package never holds the stream
    observe!(summary, batch)
end
```

Summaries from independent chunks or devices combine with [`merge!`](@ref) — the mergeable-monoid
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

**4. Two-time statistics.** For conditional and transition densities, the consumer forms lag-labelled
pairs ``\bigl(X(t), Y(t+\tau)\bigr)`` and bins them with
[`observe_pairs!`](@ref) into a
[`Histogram2D`](@ref); the snapshot buffer and lag policy
stay consumer-side. The same post-hoc rule applies — every estimand is read off the one joint summary:

```julia
joint = Histogram2D(-5.0, 5.0, 300, -5.0, 5.0, 300; lag = 1)
for (x_past, y_now) in pairs    # the consumer forms the lagged pairs at lag τ
    observe_pairs!(joint, x_past, y_now)
end

m  = conditional_mean(joint, 0.0; h = 0.3)                     # E[Y | X = 0]
s2 = conditional_variance(joint, 0.0; h = 0.3)                 # Var[Y | X = 0]
pτ = transition_density(joint, 0.0, 0.0; hx = 0.3, hy = 0.2)   # transition density p_τ(0 | 0)
```

## Where to go next

- The [Architecture](architecture.md) explains *why* this split (estimand-agnostic accumulation,
  post-hoc smoothing) is the load-bearing decision.
- The [API reference](api.md) documents every estimand — the population quantity before the estimator.
