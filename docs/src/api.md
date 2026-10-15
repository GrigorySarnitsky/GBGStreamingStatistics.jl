# API overview

The pages that follow list the complete public API, grouped by module. Every estimator is documented
**estimand first** — the population quantity, in standard probability notation, before the estimator
that approximates it from the accumulated summary.

The lower half ([`Accumulators`](api/Accumulators.md), [`Histograms`](api/Histograms.md), with the
shared smoothing kernels in [`Kernels`](api/Kernels.md)) holds the summaries and the accumulation
primitives; the upper half — [`Expectations`](api/Expectations.md), [`MarginalDensity`](api/MarginalDensity.md),
[`ConditionalExpectations`](api/ConditionalExpectations.md), [`ConditionalDensity`](api/ConditionalDensity.md)
and [`Bandwidth`](api/Bandwidth.md) — holds the estimands.
