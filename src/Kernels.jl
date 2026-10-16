# Kernels — the smoothing kernels and the kernel-moment quadrature: the cross-cutting weighting
# primitives shared by every estimator (the `density` KDE, the local-polynomial fits, the double-kernel
# conditional density, the bandwidth rules). Two families: the **symmetric** kernels (`gaussian_kernel` /
# `epanechnikov_kernel`, any callable with ``\int K = 1``; their constants — roughness ``R(K) = \int K^2``,
# second moment ``\mu_2(K) = \int u^2 K``, … — follow by quadrature `_quad`, so a user kernel needs no
# hand-coded constants) and the **asymmetric** varying-shape boundary kernels (`gamma_kernel` /
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

The **Epanechnikov** kernel ``K(u) = \\tfrac{3}{4}(1 - u^2)`` for ``|u| \\le 1`` (zero outside) is the
MSE-optimal second-order kernel (``\\int K = 1``, second moment ``\\mu_2(K) = 1/5``, roughness ``R(K) = 3/5``). Pass
it as `kernel = epanechnikov_kernel` to [`density`](@ref GBGStreamingStatistics.MarginalDensity.density);
its compact support touches only the bins within ``h`` of the query point. The
[`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) selector is **kernel-aware**, so it returns
the matched Epanechnikov scale automatically — no manual rescaling.
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


# --- asymmetric (varying-shape) kernels for bounded / positive support (Chen 1999, 2000) ---
# An asymmetric kernel is a density in the DATA argument whose shape depends on the query point x and
# bandwidth b, with support matching the data domain — so no mass leaks past the boundary and the estimate
# f̂(x) = (1/n) Σ_j n_j K(c_j; x, b) is boundary-bias-free by construction (O(b) bias uniformly, vs a
# symmetric kernel's O(h) inconsistency at an endpoint). Still a reduction over the bin counts (post-hoc,
# mergeable). The marker type routes these to the asymmetric estimator in `density`.
"""
    AsymmetricKernel(pdf, support; local_scale_weight = nothing)

A varying-shape **asymmetric kernel** for
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density) on a bounded or half-line support is
the boundary-bias-free analog of a symmetric smoothing kernel (Chen 1999, 2000). The built-ins are
[`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) for ``[0,\\infty)`` and
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) for ``[0,1]``; construct one
directly to plug in a custom boundary kernel.

# Fields

- `pdf` — the kernel, a callable `(c, x, b)` giving a density in the data argument.
- `support` — the data domain the kernel covers, a `(lower, upper)` tuple.
- `local_scale_weight` — the AMISE local-scale weight ``w(x)``, a callable, or `nothing` when unset.

# Arguments

The constructor takes the three fields (`local_scale_weight` as a keyword), the construction contract:

- `pdf(c, x, b)` — the kernel as a **density in the data argument** ``c``, whose shape is set by the
  evaluation point ``x`` and bandwidth ``b`` (normalized: ``\\int \\texttt{pdf}(c, x, b)\\,\\mathrm{d}c = 1``
  over the support, for each ``x`` and ``b``). The estimate is the binned reduction
  ``\\hat f(x) = \\tfrac{1}{n} \\sum_j n_j \\texttt{pdf}(c_j, x, b)`` over the bin counts ``n_j`` at centers
  ``c_j``, computed after accumulation from the stored counts, so estimates from merged summaries agree.
- `support` — the data domain the kernel covers, as a `(lower, upper)` tuple (e.g. `(0.0, Inf)`).
- `local_scale_weight` — the AMISE **local-scale weight** ``w(x)``, a callable `x -> w(x)`: the leading bias is
  ``\\tfrac{1}{2} b w(x) f''(x)`` and the variance ``\\propto w(x)^{-1/2}`` (``w(x) = x`` for the gamma
  kernel, ``x(1-x)`` for the beta). It is used **only** by
  [`bandwidth(…; rule = :plugin)`](@ref GBGStreamingStatistics.Bandwidth.bandwidth). Left unset
  (`nothing`, the default), `rule = :plugin` errors rather than use a wrong weight, while the ``\\hat\\sigma n^{-2/5}`` rule of thumb
(`rule = :normal`) needs no weight.

```julia
using SpecialFunctions: loggamma
# the built-in gamma kernel, as a template (Chen's modified f₂ shape); the AMISE weight
# w(x) = x serves rule = :plugin
mygamma = AsymmetricKernel(
    (c, x, b) -> c > 0 ?
        (ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b;
         exp((ρ - 1) * log(c) - c / b - ρ * log(b) - loggamma(ρ))) : 0.0,
    (0.0, Inf);
    local_scale_weight = x -> x)
# summary: a Histogram1D of the stream, accumulated with observe!
density(summary, 1.0; kernel = mygamma)
```

# References

- S. X. Chen (1999), [*Beta kernel estimators for density functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
- S. X. Chen (2000), [*Probability density function estimation using gamma kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3).
"""
struct AsymmetricKernel{F,W}
  pdf::F                            # (c, x, b) → a density in the data argument c
  support::Tuple{Float64,Float64}   # the data domain the kernel covers
  local_scale_weight::W                         # the AMISE local-scale weight w(x) (x for the gamma kernel, x(1−x)
                                    # for the beta): both the leading bias ½b·w(x)f''(x) and the variance
                                    # ∝ w(x)^{-1/2} carry it. Used by the §4.5 plug-in bandwidth.
                                    # `nothing` (no weight given) ⇒ rule=:plugin is unavailable.

  # The weight has a default — a custom kernel that only needs the rule-of-thumb bandwidth omits it — so
  # it is a keyword (RULE.API.POSITIONAL); this one constructor replaces the positional defaults.
  AsymmetricKernel(pdf, support; local_scale_weight = nothing) =
    new{typeof(pdf), typeof(local_scale_weight)}(pdf, support, local_scale_weight)
end

@inline (k::AsymmetricKernel)(c, x, b) = k.pdf(c, x, b)


# Chen's **modified** gamma kernel (f̂₂): shape ρ(x) = x/b in the interior (so the kernel mean is exactly
# x, removing the f' bias term of the naive x/b+1 form) and ρ(x) = (x/b)²/4 + 1 within 2b of the edge.
function _gamma_pdf(c, x, b)
  c > 0 || return 0.0
  ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b
  return exp((ρ - 1) * log(c) - c / b - ρ * log(b) - loggamma(ρ))
end

"""
    gamma_kernel

The **gamma kernel** (Chen 2000) estimates a density on the half-line ``[0,\\infty)`` as a varying-shape
*asymmetric* kernel: at the query point ``x`` it is the modified (bias-reduced) ``\\mathrm{Gamma}(\\rho(x), b)``
density in the data argument (Chen's ``f_2``: shape ``\\rho(x) = x/b`` in the interior — so the kernel mean
is exactly ``x``, removing the ``f'`` bias of the naive ``x/b+1`` form — and ``(x/b)^2/4 + 1`` within ``2b``
of the edge), so its support matches the data and **no mass leaks past 0**. Passed as
`kernel = gamma_kernel`, [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) returns the
boundary-bias-free estimate ``\\hat f(x) = \\tfrac{1}{n} \\sum_j n_j \\mathrm{Gamma}(c_j; \\rho(x), b)``
(``O(b)`` bias uniformly, including at the boundary where a symmetric kernel is inconsistent; the variance
self-adapts, tightening toward the edge). Use it for non-negative quantities such as droplet radii — usually
sharper than a symmetric kernel's boundary renormalization (the `support` option of `density`). The bandwidth defaults to
``\\hat\\sigma n^{-2/5}`` (the gamma-kernel rate; pass `h` to override, or the data-driven Chen plug-in
[`bandwidth(summary; rule = :plugin, kernel = gamma_kernel)`](@ref GBGStreamingStatistics.Bandwidth.bandwidth),
which lowers the ISE on skewed data). The `degree`, `adaptive` and `support` options of `density` do
not apply — the kernel is its own degree-0, boundary-free estimator.

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

The **beta kernel** (Chen 1999) estimates a density on the unit interval ``[0,1]`` as the bounded-interval
analog of [`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel): at ``x`` it is
the modified ``\\mathrm{Beta}`` density in the data argument (Chen's ``f_2``: interior shapes
``(x/b, (1-x)/b)``; the ``\\sqrt{ \\cdot }`` correction within ``2b`` of an edge replaces that shape by a
quadratic's root in the bandwidth and position), boundary-bias-free at **both** edges. Pass
`kernel = beta_kernel` with data scaled to ``[0,1]``; bandwidth default, option restrictions and the
``n^{-2/5}`` rate are as for `gamma_kernel`.

Reference: S. X. Chen (1999), [*Beta kernel estimators for density
functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
"""
const beta_kernel = AsymmetricKernel(_beta_pdf, (0.0, 1.0); local_scale_weight = x -> x * (1 - x))


# The kernel-weighting atom — the one place the two kernel protocols live. Given a bin/data location `c`,
# a query point `x`, and a bandwidth `h`, return the per-location smoothing weight, **normalized** (the
# 1/h is folded in for a symmetric kernel, so every downstream readout — a count sum, a Nadaraya–Watson
# ratio, the local-polynomial jet — is uniform, and the symmetric/asymmetric split never re-appears).
# A symmetric kernel is the scaled difference (1/h)K((c−x)/h); an asymmetric kernel is already a normalized
# density in its data argument, K(c; x, h) (no 1/h). Internal (reached by explicit `using`).
@inline _kernel_weight(c, x, kernel, h) = kernel((c - x) / h) / h
@inline _kernel_weight(c, x, kernel::AsymmetricKernel, h) = kernel(c, x, h)


end # module Kernels
