# ConditionalDensity

The conditional density [`conditional_density`](@ref) (``f_{Y\mid X}``) and the two-time transition density [`transition_density`](@ref) (``p_\tau(y\mid x)``), by the double-kernel estimator over a separate conditioning `kernelx` and response `kernely` (the asymmetric bounded-response kernels among them; see the [Kernels](Kernels.md) page) — each with its pointwise standard error and a `_with_standard_error` reader returning estimate and error together.

```@docs
conditional_density
conditional_density_standard_error
conditional_density_with_standard_error
transition_density
transition_density_standard_error
transition_density_with_standard_error
```
