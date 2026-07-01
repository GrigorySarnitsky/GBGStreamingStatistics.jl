# Expectations

The power-sum quadrature readouts (§4.B) — a smooth functional of ``X`` read straight off the stored
power-sums, with no kernel and no query point. `expectation` estimates the expectation
``\mathbb{E}[g(X)]`` of a general smooth `g`; `expectation_standard_error` its sampling standard error;
and `pushforward_density` the density of the transformed variable ``g(X)``. All three route through the
one shared contraction atom that the recovered moments (`mean` / `var` / `central_moment`) also use, so
those are the ``g = x^p`` case; a general `g`'s per-bin coefficients come from a dependency-free local
polynomial interpolant (exact on monomials), and the pushforward is `expectation` of the composed kernel.

```@docs
expectation
expectation_standard_error
pushforward_density
```

## Recovered moments

The global marginal moments — the ``g = x^p`` case — recovered post hoc from the power-sums by the shift
formula (exact), extending `Statistics`.

```@docs
mean(::Histogram1D)
var(::Histogram1D)
std(::Histogram1D)
central_moment
```
