# MarginalDensity

The marginal density ``\equiv f_X``: the flagship [`density`](@ref), its pointwise variability band
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

The moment-fitted parametric models for the Hjort–Glad semiparametric estimator `density(…; start = …)`,
which multiplies the start by a KDE of its correction. The built-ins are `normal_start`, `lognormal_start`,
`gamma_start`, `exponential_start`, and `inverse_gaussian_start`; a custom model is a
[`ParametricStart`](@ref) — supply the moment fit, the log-density, and the support.

```@docs
normal_start
lognormal_start
gamma_start
exponential_start
inverse_gaussian_start
ParametricStart
```
