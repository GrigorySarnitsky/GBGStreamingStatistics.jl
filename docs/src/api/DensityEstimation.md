# DensityEstimation

The marginal density ``\equiv f_X`` (§4.1): the flagship [`density`](@ref) and its pointwise variability
band, the smoothing kernels it accepts, and the parametric starts for its semiparametric mode.

## Density

```@docs
density
density_standard_error
```

## Kernels

The smoothing kernel passed as `density(…; kernel = …)` — the symmetric `gaussian_kernel` /
`epanechnikov_kernel`, and the bounded/positive-support **asymmetric** `gamma_kernel` / `beta_kernel`,
which are boundary-bias-free for data on ``[0,\infty)`` or ``[0,1]``.

```@docs
gaussian_kernel
epanechnikov_kernel
gamma_kernel
beta_kernel
```

## Parametric starts

The moment-fitted parametric models for the Hjort–Glad semiparametric estimator `density(…; start = …)`,
which multiplies the start by a KDE of its correction.

```@docs
normal_start
lognormal_start
gamma_start
exponential_start
inverse_gaussian_start
```
