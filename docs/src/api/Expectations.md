# Expectations

The functions on this page read a smooth functional of ``X`` off the power-sums a [`Histogram1D`](@ref)
of its samples stores (the summed powers of each bin's sample offsets from its center), with no kernel
and no query point. The function [`expectation`](@ref) estimates the expectation
``\mathbb{E}[g(X)]`` of a smooth `g`; [`expectation_standard_error`](@ref) its sampling standard error
(paired with the estimate by [`expectation_with_standard_error`](@ref)); and [`pushforward_density`](@ref) the
density of the transformed variable ``g(X)``.

All three use the same weighted sum of the power-sums as the sample moments (`mean` / `var` / `std` /
[`central_moment`](@ref)), so those are the ``g = x^p`` case; a general `g`'s per-bin coefficients come
from a local polynomial interpolant (exact on monomials), and the pushforward is the `expectation` of the
kernel composed with `g`.

```@docs
expectation
expectation_standard_error
expectation_with_standard_error
pushforward_density
pushforward_density_standard_error
pushforward_density_with_standard_error
```

## Recovered moments

The global marginal moments — the ``g = x^p`` case — are recovered exactly from the power-sums by the
shift formula, the binomial re-expansion of each bin's powers about its center, and extend `Statistics`.

```@docs
mean(::Histogram1D)
var(::Histogram1D)
std(::Histogram1D)
central_moment
```
