# Expectations

The power-sum quadrature readouts read a smooth functional of ``X`` straight off the stored
power-sums, with no kernel and no query point. The reader [`expectation`](@ref) estimates the expectation
``\mathbb{E}[g(X)]`` of a general smooth `g`; [`expectation_standard_error`](@ref) its sampling standard error
(paired with the estimate in one call by [`expectation_with_standard_error`](@ref)); and [`pushforward_density`](@ref) the
density of the transformed variable ``g(X)``. All three are the same power-sum contraction the recovered
moments (`mean` / `var` / `std` / `central_moment`) use, so those are the ``g = x^p`` case; a general
`g`'s per-bin coefficients come from a dependency-free local
polynomial interpolant (exact on monomials), and the pushforward is `expectation` of the composed kernel.

```@docs
expectation
expectation_standard_error
expectation_with_standard_error
pushforward_density
pushforward_density_standard_error
pushforward_density_with_standard_error
```

## Recovered moments

The global marginal moments — the ``g = x^p`` case — are recovered post hoc from the power-sums by the
shift formula, exactly, extending `Statistics`.

```@docs
mean(::Histogram1D)
var(::Histogram1D)
std(::Histogram1D)
central_moment
```
