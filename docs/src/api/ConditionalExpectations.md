# ConditionalExpectations

The conditional expectation ``\mathbb{E}[g(Y) \mid X = x]`` of a smooth functional of the response
([`conditional_expectation`](@ref)) with its sampling standard error
([`conditional_expectation_standard_error`](@ref), paired with the estimate by
[`conditional_expectation_with_standard_error`](@ref)), the conditional mean [`conditional_mean`](@ref)
(``m(x)``) and variance [`conditional_variance`](@ref) (``\sigma^2(x)``, with its residual-based
standard error) as its ``g = \mathrm{identity}`` / ``y^2`` cases, and the conditional pushforward density
[`conditional_pushforward_density`](@ref) (``f_{g(Y) \mid X}``) of a transformed response. All are
local-polynomial regressions in ``x`` of a ``y``-quadrature response (the per-column expectation of the
functional, read off the stored response power-sums).

```@docs
conditional_expectation
conditional_expectation_standard_error
conditional_expectation_with_standard_error
conditional_mean
conditional_variance
conditional_variance_standard_error
conditional_variance_with_standard_error
conditional_pushforward_density
conditional_pushforward_density_standard_error
conditional_pushforward_density_with_standard_error
```
