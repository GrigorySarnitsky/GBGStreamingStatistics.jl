# Internals

Shared building blocks *below* the public estimator API — documented here for development, review, and
**extension** (building a new local-polynomial estimator on the same primitive without reimplementing the
solve — the architecture's shared core, [Architecture](architecture.md) §11).

!!! warning "Not the stable public API"
    These bindings are **not exported** and **not part of the stable surface** — reach them by the
    qualified name (`GBGStreamingStatistics.LocalPolynomial.local_fit`) or an explicit
    `using ..LocalPolynomial: local_fit`, and expect signatures to change. Routine analysis uses the
    exported estimands ([`density`](@ref GBGStreamingStatistics.DensityEstimation.density),
    [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean),
    [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth), …), which take the `Histogram` summary
    types; `local_fit` works one level lower, on raw recovered arrays.

## The local-polynomial jet — `LocalPolynomial.local_fit`

The single linear primitive on which the local smoothers are built (the formalization:
[Formalization](formalization.md), "Local-polynomial jet primitive"). `conditional_mean` is its
``\beta_0``; the local-likelihood `density` reuses its solve as the iteratively reweighted inner step.

```@docs
GBGStreamingStatistics.LocalPolynomial.local_fit
```
