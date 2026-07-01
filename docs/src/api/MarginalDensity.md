# MarginalDensity

The marginal density ``\equiv f_X``: the [`density`](@ref) estimator, its pointwise variability band
[`density_standard_error`](@ref) and the paired [`density_with_standard_error`](@ref), and the parametric
starts for its semiparametric mode. `density` accepts any symmetric or bounded-support kernel from the
[Kernels](Kernels.md) page.

## Density

```@docs
density
density_standard_error
density_with_standard_error
```

## Parametric starts

The moment-fitted [`ParametricStart`](@ref) models for the Hjort–Glad semiparametric
`density(…; start = …)`, each documented below.

```@docs
normal_start
lognormal_start
gamma_start
exponential_start
inverse_gaussian_start
ParametricStart
```
