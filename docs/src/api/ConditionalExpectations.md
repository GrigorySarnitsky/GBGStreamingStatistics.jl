# ConditionalExpectations

The conditional expectation ``\mathbb{E}[g(Y) \mid x]`` of a smooth functional of the response
(`conditional_expectation`) with its sampling standard error (`conditional_expectation_standard_error`,
paired with the estimate in one call by `conditional_expectation_with_standard_error`), the conditional
mean `conditional_mean` ``\equiv m(x)`` and variance `conditional_variance`
``\equiv \sigma^2(x)`` as its ``g = \mathrm{identity}`` / ``y^2`` cases, and the conditional pushforward
density `conditional_pushforward_density` ``\equiv f_{g(Y) \mid X}`` of a transformed response (§4.2). All
by local-polynomial regression in ``x`` of a per-column ``y``-quadrature — composing both atoms
(kernel-weighting in ``x`` ∘ power-sum quadrature in ``y``). The estimand-first reference module.

```@docs
conditional_expectation
conditional_expectation_standard_error
conditional_expectation_with_standard_error
conditional_mean
conditional_variance
conditional_pushforward_density
```
