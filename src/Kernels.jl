# Kernels — the symmetric smoothing kernels and the kernel-moment quadrature: the cross-cutting
# weighting primitives shared by every estimator (the `density` KDE, the local-polynomial fits, the
# double-kernel conditional density, the bandwidth rules). A kernel is any callable with ``\int K = 1``;
# its constants (roughness ``R(K) = \int K^2``, second moment ``\mu_2(K) = \int u^2 K``, …) follow by
# quadrature (`_quad`), so a user-supplied kernel needs no hand-coded constants.

module Kernels

export gaussian_kernel, epanechnikov_kernel


"""
    gaussian_kernel(u)

The standard Gaussian smoothing kernel ``K(u) = e^{-u^2/2}/\\sqrt{2\\pi}`` (``\\int K = 1``,
``\\mu_2(K) = 1``, roughness ``R(K) = 1/(2\\sqrt{\\pi})``) — the
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density) default.
"""
@inline gaussian_kernel(u) = exp(-u^2 / 2) / sqrt(2π)


"""
    epanechnikov_kernel(u)

The **Epanechnikov** kernel ``K(u) = \\tfrac{3}{4}(1 - u^2)`` for ``|u| \\le 1`` (zero outside) — the
MSE-optimal second-order kernel (``\\int K = 1``, ``\\mu_2(K) = 1/5``, roughness ``R(K) = 3/5``). Pass
it as `kernel = epanechnikov_kernel`; its compact support touches only the bins within ``h`` of the
query point. [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) is **kernel-aware**, so it
returns the matched Epanechnikov ``h \\approx 2.34\\,\\hat\\sigma\\,n^{-1/5}`` (``2.21\\times`` the
Gaussian ``h``) automatically — no manual rescaling.
"""
@inline epanechnikov_kernel(u) = abs(u) < 1 ? 0.75 * (1 - u^2) : 0.0


# Kernel-moment quadrature: integrate `g` over [−L, L] by the trapezoidal rule on a fixed fine grid, so
# any user-supplied kernel's constants (``R(K)``, ``\mu_2(K)``, ``\mu_4(K)``, the CDF, …) follow by
# quadrature. The shipped kernels sit well inside the range (Gaussian ≈ ±8, Epanechnikov ±1); this runs
# once per `density` / `bandwidth` call, off the hot path. Internal (reached by explicit `using`).
function _quad(g, L::Float64 = 16.0, npts::Int = 32_769)
  Δ = 2L / (npts - 1)
  acc = 0.5 * (g(-L) + g(L))
  @inbounds for i in 1:npts-2
    acc += g(-L + i * Δ)
  end
  return acc * Δ
end


end # module Kernels
