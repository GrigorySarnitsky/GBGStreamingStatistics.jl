# MarginalDensity

The marginal density ``\equiv f_X`` (§4.1): the flagship [`density`](@ref) and its pointwise variability
band, the smoothing kernels it accepts, and the parametric starts for its semiparametric mode.

## Density

```@docs
density
density_standard_error
```

## Asymmetric kernels

The smoothing kernel passed as `density(…; kernel = …)`. The **symmetric** kernels (`gaussian_kernel` /
`epanechnikov_kernel`, or any callable with ``\int K = 1``) live on the Kernels page. The
bounded/positive-support **asymmetric** kernels are [`AsymmetricKernel`](@ref) values — `gamma_kernel`
for ``[0,\infty)``, `beta_kernel` for ``[0,1]``, or a custom one — boundary-bias-free by construction.

```@docs
gamma_kernel
beta_kernel
AsymmetricKernel
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
