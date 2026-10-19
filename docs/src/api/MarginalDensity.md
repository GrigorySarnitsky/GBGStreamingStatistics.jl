# MarginalDensity

The marginal density ``\equiv f_X`` is estimated by [`density`](@ref); its pointwise standard error is
[`density_standard_error`](@ref), [`density_with_standard_error`](@ref) returns both, and the
parametric starts serve its semiparametric mode. The estimator `density` accepts any symmetric or bounded-support kernel from the
[Kernels](Kernels.md) page.

## Density

```@docs
density
density_standard_error
density_with_standard_error
```

## Parametric starts

The moment-fitted [`ParametricStart`](@ref) models serve the Hjort–Glad semiparametric
`density(…; start = …)`; each is documented below.

```@docs
normal_start
lognormal_start
gamma_start
exponential_start
inverse_gaussian_start
ParametricStart
```
