# Kernels

The smoothing kernels shared across the estimators, in two families. The **symmetric** kernels are any
callable with ``\int K = 1`` (the kernel-aware [`bandwidth`](@ref) computes the roughness ``R(K)`` and
second moment ``\mu_2(K)`` by quadrature, so a custom symmetric kernel needs no manual constants). The
**asymmetric** varying-shape kernels — `gamma_kernel` for ``[0,\infty)``, `beta_kernel` for ``[0,1]``, or
a custom [`AsymmetricKernel`](@ref) — match a bounded / positive support and are boundary-bias-free by
construction (Chen 1999, 2000). Both families are passed to [`density`](@ref) as `kernel = …`.

## Symmetric

```@docs
gaussian_kernel
epanechnikov_kernel
```

## Asymmetric (bounded / positive support)

```@docs
gamma_kernel
beta_kernel
AsymmetricKernel
```
