# API overview

These pages list the public API by module. The package reads a stream of samples once, folding each
batch into a binned summary and discarding it. Every estimator is documented estimand first — the
population quantity, in standard probability notation, before the estimator that approximates it from
that summary.

The accumulation modules ([`Accumulators`](api/Accumulators.md), [`Histograms`](api/Histograms.md), with the
shared smoothing kernels in [`Kernels`](api/Kernels.md)) hold the summaries and the accumulation
primitives; the estimator modules — [`Expectations`](api/Expectations.md), [`MarginalDensity`](api/MarginalDensity.md),
[`ConditionalExpectations`](api/ConditionalExpectations.md), [`ConditionalDensity`](api/ConditionalDensity.md)
and [`Bandwidth`](api/Bandwidth.md) — hold the estimators.
