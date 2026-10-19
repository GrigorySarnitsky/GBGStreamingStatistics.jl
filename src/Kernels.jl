# Kernels — the smoothing kernels and the kernel-moment quadrature: the cross-cutting weighting
# primitives shared by every estimator (the `density` KDE, the local-polynomial fits, the double-kernel
# conditional density, the bandwidth rules). Two families: the symmetric kernels (`gaussian_kernel` /
# `epanechnikov_kernel`, any callable with ∫K = 1; their constants — roughness R(K) = ∫K²,
# second moment μ₂(K) = ∫u²K, … — follow by quadrature `_quad`, so a user kernel needs no
# hand-coded constants) and the asymmetric varying-shape boundary kernels (`gamma_kernel` /
# `beta_kernel`, the `AsymmetricKernel` type) for bounded / positive support (Chen 1999, 2000).

module Kernels

using SpecialFunctions: loggamma, logbeta

export gaussian_kernel, epanechnikov_kernel, AsymmetricKernel, gamma_kernel, beta_kernel


"""
    gaussian_kernel(u)

The standard Gaussian smoothing kernel ``K(u) = e^{-u^2/2}/\\sqrt{2\\pi}`` (``\\int K = 1``,
second moment ``\\mu_2(K) = 1``, roughness ``R(K) = 1/(2\\sqrt{\\pi})``) is the
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density) default.
"""
@inline gaussian_kernel(u) = exp(-u^2 / 2) / sqrt(2π)


"""
    epanechnikov_kernel(u)

The Epanechnikov kernel ``K(u) = \\tfrac{3}{4}(1 - u^2)`` for ``|u| \\le 1`` (zero outside) is the
MSE-optimal second-order kernel (``\\int K = 1``, second moment ``\\mu_2(K) = 1/5``, roughness ``R(K) = 3/5``). Pass
it as `kernel = epanechnikov_kernel` to [`density`](@ref GBGStreamingStatistics.MarginalDensity.density);
its compact support touches only the bins within the bandwidth ``h`` of the query point. The
[`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) selector is kernel-aware, so it returns
the matched Epanechnikov scale directly.
"""
@inline epanechnikov_kernel(u) = abs(u) < 1 ? 0.75 * (1 - u^2) : 0.0


# Kernel-moment quadrature: integrate `g` over `[-L, L]` (`L` the half-width) by the trapezoidal rule on a
# fixed fine grid, so any user-supplied kernel's constants (R(K), μ₂(K), μ₄(K), the CDF, …) follow by
# quadrature. The shipped kernels sit well inside the range (Gaussian ≈ ±8, Epanechnikov ±1); this runs
# once per `density` / `bandwidth` call, off the hot path. Internal (reached by explicit `using`).
function _quad(g; L::Float64 = 16.0, npts::Int = 32_769)
  Δ = 2L / (npts - 1)
  acc = 0.5 * (g(-L) + g(L))
  @inbounds for i in 1:npts-2
    acc += g(-L + i * Δ)
  end
  return acc * Δ
end


# --- asymmetric (varying-shape) kernels for bounded / positive support (Chen 1999, 2000) ---
# An asymmetric kernel is a density in the DATA argument whose shape depends on the query point x and
# bandwidth b, with support matching the data domain — so no mass leaks past the boundary and the estimate
# f̂_X(x) = (1/n) Σ_j n_j K(c_j; x, b) is boundary-bias-free by construction (O(b) bias uniformly, vs a
# symmetric kernel's O(1) bias at an endpoint). Still a reduction over the bin counts (post-hoc,
# mergeable). The marker type routes these to the asymmetric estimator in `density`.
"""
    AsymmetricKernel(pdf, support; local_scale_weight = nothing)

A varying-shape *asymmetric kernel* for
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density) on a bounded or half-line support is
the boundary-bias-free analog of a symmetric smoothing kernel (Chen 1999, 2000). The built-ins are
[`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) for ``[0,\\infty)`` and
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) for ``[0,1]``; construct one
directly to plug in a custom boundary kernel and its local-scale weight ``w(x)`` for the plug-in
bandwidth.

# Arguments

The constructor takes the three fields, the last as a keyword:

- `pdf(v, x, b)` — the kernel as a density in the data argument ``v``, whose shape is set by the
  evaluation point ``x`` and bandwidth ``b`` (normalized: ``\\int \\texttt{pdf}(v, x, b)\\,\\mathrm{d}v = 1``
  over the support, for each ``x`` and ``b``). `density` sums it over the bin counts at the bin
  locations, so merged summaries give the same estimate.
- `support` — the data domain the kernel covers, as a `(lower, upper)` tuple (e.g. `(0.0, Inf)`).
- `local_scale_weight` — the AMISE local-scale weight ``w(x)``, a callable `x -> w(x)`: the leading bias is
  ``\\tfrac{1}{2} b w(x) f_X''(x)`` and the interior variance ``\\propto w(x)^{-1/2}`` (``w(x) = x`` for the gamma
  kernel, ``x(1-x)`` for the beta). Only
  [`bandwidth(…; rule = :plugin)`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) uses it. If it is
  unset (`nothing`, the default), `rule = :plugin` errors, while the
  ``\\hat\\sigma_X n^{-2/5}`` rule of thumb (`rule = :normal`, from the sample standard deviation) needs no weight.

```julia
using SpecialFunctions: loggamma
# the built-in gamma kernel, as a template (the shape of Chen's bias-reduced
# estimator); the AMISE weight w(x) = x serves rule = :plugin
mygamma = AsymmetricKernel(
    (v, x, b) -> v > 0 ?
        (ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b;
         exp((ρ - 1) * log(v) - v / b - ρ * log(b) - loggamma(ρ))) : 0.0,
    (0.0, Inf);
    local_scale_weight = x -> x)
# summary: a Histogram1D of the stream, accumulated with observe!
density(summary, 1.0; kernel = mygamma)
```

# References

- S. X. Chen (1999), [*Beta kernel estimators for density functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
- S. X. Chen (2000), [*Probability density function estimation using gamma kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3).

# Fields

- `pdf` — the kernel, a callable `(v, x, b)` giving a density in the data argument.
- `support` — the data domain the kernel covers, a `(lower, upper)` tuple.
- `local_scale_weight` — the AMISE local-scale weight, a callable, or `nothing` when unset.
"""
struct AsymmetricKernel{F,W}
  pdf::F                            # (v, x, b) → a density in the data argument v
  support::Tuple{Float64,Float64}   # the data domain the kernel covers
  local_scale_weight::W                         # the AMISE local-scale weight w(x) (x for the gamma kernel, x(1−x)
                                    # for the beta): both the leading bias ½b·w(x)f_X''(x) and the variance
                                    # ∝ w(x)^{-1/2} carry it. Used by the §4.5 plug-in bandwidth.
                                    # `nothing` (no weight given) ⇒ rule=:plugin is unavailable.

  # The weight has a default — a custom kernel that only needs the rule-of-thumb bandwidth omits it — so
  # it is a keyword (RULE.API.POSITIONAL); this one constructor replaces the positional defaults.
  AsymmetricKernel(pdf, support; local_scale_weight = nothing) =
    new{typeof(pdf), typeof(local_scale_weight)}(pdf, support, local_scale_weight)
end

@inline (k::AsymmetricKernel)(c, x, b) = k.pdf(c, x, b)


# Chen's modified gamma kernel (the kernel of his bias-reduced estimator): shape ρ(x) = x/b in the interior (so the kernel
# mean is exactly x, removing the f_X' bias term of the naive x/b+1 form) and ρ(x) = (x/b)²/4 + 1 within 2b of the edge.
function _gamma_pdf(c, x, b)
  c > 0 || return 0.0
  ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b
  return exp((ρ - 1) * log(c) - c / b - ρ * log(b) - loggamma(ρ))
end

"""
    gamma_kernel

The gamma kernel (Chen 2000) estimates a density on the half-line ``[0,\\infty)`` as an asymmetric
kernel: at the query point ``x`` it is the modified ``\\mathrm{Gamma}(\\rho(x), b)`` density in the
data argument, with the bandwidth (the keyword `h`) as its scale. Chen's bias-reduced shape is
``\\rho(x) = x/b`` in the interior, so the kernel mean is exactly ``x``, removing the slope ``f_X'``
from the bias of the naive ``x/b+1`` shape, and ``(x/b)^2/4 + 1`` within ``2b`` of the edge.

With `kernel = gamma_kernel`, [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) returns
the estimate from the bin counts at the bin locations, with ``O(b)`` bias uniformly, also at the boundary, where a symmetric kernel is
inconsistent; the variance grows toward the edge. Use it for non-negative data, usually sharper than
the `support` renormalization of `density`. The bandwidth defaults to ``\\hat\\sigma_X n^{-2/5}``, the
sample standard deviation times the sample-size rate; override it with `h`, or use the Chen plug-in
[`bandwidth(summary; rule = :plugin, kernel = gamma_kernel)`](@ref GBGStreamingStatistics.Bandwidth.bandwidth),
which lowers the ISE on skewed data. The local-polynomial `degree`, variable-bandwidth `adaptive` and
`support` options of `density` do not apply: the kernel is its own degree-0 estimator.

Reference: S. X. Chen (2000), [*Probability density function estimation using gamma
kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3).
"""
const gamma_kernel = AsymmetricKernel(_gamma_pdf, (0.0, Inf); local_scale_weight = x -> x)


# Chen's modified beta shape: interior (x/b, (1-x)/b); a √-correction within 2b of an edge (mean ≈ x).
_beta_rho(x, b) = 2b^2 + 2.5 - sqrt(4b^4 + 6b^2 + 2.25 - x^2 - x / b)

function _beta_pdf(c, x, b)
  (0 < c < 1) || return 0.0
  p = x < 2b ? _beta_rho(x, b) : x / b
  q = x > 1 - 2b ? _beta_rho(1 - x, b) : (1 - x) / b
  return exp((p - 1) * log(c) + (q - 1) * log1p(-c) - logbeta(p, q))
end

"""
    beta_kernel

The beta kernel (Chen 1999) estimates a density on ``[0,1]`` as the bounded-interval
analog of [`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel): at ``x`` it is
the modified ``\\mathrm{Beta}`` density in the data argument, boundary-bias-free at both edges. Its
interior shapes are ``(x/b, (1-x)/b)``, over the bandwidth `h`; a ``\\sqrt{ \\cdot }``
correction within ``2b`` of an edge replaces that side's shape by a quadratic's root in bandwidth
and edge distance. Pass `kernel = beta_kernel` with data scaled to ``[0,1]``;
bandwidth default, option restrictions and the ``n^{-2/5}`` rate are as for `gamma_kernel`.

Reference: S. X. Chen (1999), [*Beta kernel estimators for density
functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
"""
const beta_kernel = AsymmetricKernel(_beta_pdf, (0.0, 1.0); local_scale_weight = x -> x * (1 - x))


# The kernel-weighting atom — the one place the two kernel protocols live. Given a bin/data location `c`,
# a query point `x`, and a bandwidth `h`, return the per-location smoothing weight, normalized (the
# 1/h is folded in for a symmetric kernel, so every downstream readout — a count sum, a Nadaraya–Watson
# ratio, the local-polynomial jet — is uniform, and the symmetric/asymmetric split never re-appears).
# A symmetric kernel is the scaled difference (1/h)K((c−x)/h); an asymmetric kernel is already a normalized
# density in its data argument, K(c; x, b), the bandwidth h serving as its scale b (no 1/h). Internal (reached by explicit `using`).
@inline _kernel_weight(c, x, kernel, h) = kernel((c - x) / h) / h
@inline _kernel_weight(c, x, kernel::AsymmetricKernel, h) = kernel(c, x, h)


end # module Kernels
