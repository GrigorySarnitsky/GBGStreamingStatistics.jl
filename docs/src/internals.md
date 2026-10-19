# Internals

These are the shared building blocks *below* the public estimator API, documented here for development,
review, and extension (a new local-polynomial estimator on the local-polynomial fit `local_fit`, a kernel-weighted
least-squares solve, without reimplementing it; see the [Architecture](architecture.md)).

!!! warning "Not the stable public API"
    These bindings are not exported and not part of the stable surface — reach them by the
    qualified name (`GBGStreamingStatistics.LocalPolynomial.local_fit`) or an explicit
    `using ..LocalPolynomial: local_fit`, and expect signatures to change. Routine analysis uses the
    exported estimands ([`density`](@ref), [`conditional_mean`](@ref), [`bandwidth`](@ref), …), which take a [`Histogram1D`](@ref) /
    [`Histogram2D`](@ref) summary (the binned record of the stream); `local_fit` works one level
    lower, on arrays read back from one.

## The local-polynomial jet — `LocalPolynomial.local_fit`

This is the single linear primitive on which the local smoothers are built (the formalization:
[Formalization](formalization.md), "Local-polynomial jet primitive"). The estimator `conditional_mean` is its
fitted value ``\hat\beta_0``, the leading coefficient of the jet (the fitted local Taylor coefficients); the local-likelihood `density` reuses its solve as the iteratively reweighted inner step.

```@docs
GBGStreamingStatistics.LocalPolynomial.local_fit
```
