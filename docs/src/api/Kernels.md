# Kernels

The smoothing kernels shared across the estimators come in two families; every estimator applies one to a
`summary`, the binned histogram of the stream. The symmetric kernels are any callable with ``\int K = 1`` (the kernel-aware [`bandwidth`](@ref) computes ``R(K)`` and ``\mu_2(K)`` by
quadrature, so a custom symmetric kernel needs no manual constants). The
asymmetric varying-shape kernels — [`gamma_kernel`](@ref) for ``[0,\infty)``, [`beta_kernel`](@ref) for ``[0,1]``, or
a custom [`AsymmetricKernel`](@ref) — match a bounded / positive support and are boundary-bias-free by
construction (Chen 1999, 2000). Both families feed the estimators as `kernel = …`: the marginal [`density`](@ref), the
[`conditional_density`](@ref) (as `kernely`), the [`conditional_expectation`](@ref), and the
`bandwidth` rules.

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
