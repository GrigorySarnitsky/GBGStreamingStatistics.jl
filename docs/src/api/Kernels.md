# Kernels

The symmetric smoothing kernels shared across the estimators — any callable with ``\int K = 1`` works
(the kernel-aware [`bandwidth`](@ref) computes its roughness ``R(K)`` and second moment ``\mu_2(K)`` by
quadrature, so a custom symmetric kernel needs no manual constants). Two ship; the bounded/positive-support
**asymmetric** kernels live with `MarginalDensity`.

```@docs
gaussian_kernel
epanechnikov_kernel
```
