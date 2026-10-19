# ConditionalDensity

The conditional density [`conditional_density`](@ref) (``f_{Y\mid X}``) and the two-time transition density [`transition_density`](@ref) (``p_\tau(y\mid x)``, its subscript the lag) are estimated by a double kernel, a separate conditioning `kernel_x` and response `kernel_y` (see the [Kernels](Kernels.md) page). Each has its pointwise standard error and a `_with_standard_error` reader returning both.

```@docs
conditional_density
conditional_density_standard_error
conditional_density_with_standard_error
transition_density
transition_density_standard_error
transition_density_with_standard_error
```
