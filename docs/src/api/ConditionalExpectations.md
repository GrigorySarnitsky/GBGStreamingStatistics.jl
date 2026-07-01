# ConditionalExpectations

The conditional expectation ``\mathbb{E}[g(Y) \mid x]`` of a smooth functional of the response
([`conditional_expectation`](@ref)) with its sampling standard error
([`conditional_expectation_standard_error`](@ref), paired with the estimate in one call by
[`conditional_expectation_with_standard_error`](@ref)), the conditional mean [`conditional_mean`](@ref)
``\equiv m(x)`` and variance [`conditional_variance`](@ref) ``\equiv \sigma^2(x)`` as its
``g = \mathrm{identity}`` / ``y^2`` cases, and the conditional pushforward density
[`conditional_pushforward_density`](@ref) ``\equiv f_{g(Y) \mid X}`` of a transformed response. All
by local-polynomial regression in ``x`` of a ``y``-quadrature response.

```@docs
conditional_expectation
conditional_expectation_standard_error
conditional_expectation_with_standard_error
conditional_mean
conditional_variance
conditional_pushforward_density
```
