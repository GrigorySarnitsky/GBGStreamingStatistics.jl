# API overview

The complete public API, grouped by module on the pages that follow. Every estimator is documented
**estimand first** — the population quantity, in standard probability notation, before the estimator
that approximates it from the accumulated summary.

The lower half ([`Accumulators`](api/Accumulators.md), [`Histograms`](api/Histograms.md), with the
shared smoothing kernels in [`Kernels`](api/Kernels.md)) holds the summaries and the accumulation
primitives; the upper half ([`Expectations`](api/Expectations.md) through
[`Bandwidth`](api/Bandwidth.md)) holds the estimands.
