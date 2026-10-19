# ConditionalExpectations

The module estimates the conditional expectation ``\mathbb{E}[g(Y) \mid X = x]`` of a smooth functional of
the response ([`conditional_expectation`](@ref)) with its standard error
([`conditional_expectation_standard_error`](@ref), both at once from
[`conditional_expectation_with_standard_error`](@ref)); the conditional mean [`conditional_mean`](@ref)
(``m(x)``) and variance [`conditional_variance`](@ref) (``\sigma^2(x)``, with standard error), the
mean its ``g = \mathrm{identity}`` case and the ``y^2`` case the variance's plug-in ingredient; and the
conditional pushforward density [`conditional_pushforward_density`](@ref) (``f_{g(Y) \mid X}``). All are
local-polynomial regressions in ``x`` of a ``y``-quadrature response: the functional's mean per
conditioning bin of the joint histogram, from its stored sums of the response and its square.

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
