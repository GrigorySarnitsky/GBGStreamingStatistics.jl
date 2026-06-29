# DensityEstimation — §4.1: the marginal density f_X. A single `density` spans the local-log-polynomial
# degree (0 = KDE, ≥1 = local likelihood) crossed with a fixed or Abramson-variable bandwidth, plus
# pluggable kernels (the symmetric gaussian/epanechnikov and the bounded-support asymmetric gamma/beta),
# the Jones boundary `bounds` renormalization, the Hjort–Glad parametric `start`, exact `normalize`, and
# the `density_standard_error` band. Two stages: set the bandwidth, then fit. The bandwidth is
# *sample-point* (per-bin h_j) for the KDE and *balloon* (one h(x) at the query) for the local-log fit.
# Written estimand-first (`docs/ARCHITECTURE.md` §2, §4.1).

module DensityEstimation

using LinearAlgebra
using SpecialFunctions: loggamma, logbeta

using ..Histograms: Histogram1D, _nbins, _bincenter, _total, _sample_mean, _sample_variance, _sample_std

export density, density_standard_error, gaussian_kernel, epanechnikov_kernel,
       AsymmetricKernel, gamma_kernel, beta_kernel,
       ParametricStart, normal_start, lognormal_start, gamma_start, exponential_start, inverse_gaussian_start


"""
    gaussian_kernel(u)

The standard Gaussian smoothing kernel ``K(u) = e^{-u^2/2}/\\sqrt{2\\pi}`` (``\\int K = 1``,
``\\mu_2(K) = 1``, roughness ``R(K) = 1/(2\\sqrt{\\pi})``) — the
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density) default.
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


# --- asymmetric (varying-shape) kernels for bounded / positive support (Chen 1999, 2000) ---
# An asymmetric kernel is a density in the DATA argument whose shape depends on the query point x and
# bandwidth b, with support matching the data domain — so no mass leaks past the boundary and the estimate
# f̂(x) = (1/n) Σ_j n_j K(c_j; x, b) is boundary-bias-free by construction (O(b) bias uniformly, vs a
# symmetric kernel's O(h) inconsistency at an endpoint). Still a reduction over the bin counts (post-hoc,
# mergeable). The marker type routes these to the asymmetric estimator in `density`.
"""
    AsymmetricKernel(pdf, support)
    AsymmetricKernel(pdf, support, weight)

A varying-shape **asymmetric kernel** for
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density) on a bounded or half-line support —
the boundary-bias-free analog of a symmetric smoothing kernel (Chen 1999, 2000). The built-ins are
[`gamma_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gamma_kernel) for ``[0,\\infty)`` and
[`beta_kernel`](@ref GBGStreamingStatistics.DensityEstimation.beta_kernel) for ``[0,1]``; construct one
directly to plug in a custom boundary kernel. The three fields are the construction contract:

- `pdf(c, x, b)` — the kernel as a **density in the data argument** ``c``, whose shape is set by the
  evaluation point ``x`` and bandwidth ``b`` (normalized: ``\\int \\texttt{pdf}(c, x, b)\\,\\mathrm{d}c = 1``
  over the support, for each ``x`` and ``b``). The estimate is the binned reduction
  ``\\hat f(x) = \\tfrac1n \\sum_j n_j\\,\\texttt{pdf}(c_j, x, b)`` over the bin counts ``n_j`` at centers
  ``c_j`` — post hoc and mergeable.
- `support` — the data domain the kernel covers, as a `(lower, upper)` tuple (e.g. `(0.0, Inf)`).
- `weight(x)` — the AMISE **local-scale weight** ``w(x)``: the leading bias is
  ``\\tfrac12 b\\,w(x)\\,f''(x)`` and the variance ``\\propto w(x)^{-1/2}`` (``w(x) = x`` for the gamma
  kernel, ``x(1-x)`` for the beta). It is used **only** by
  [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth)`(…; rule = :plugin)`. The two-argument
  constructor leaves it unset (`nothing`): `rule = :plugin` then errors rather than silently use a wrong
  weight, while the ``\\hat\\sigma\\,n^{-2/5}`` rule of thumb (`rule = :normal`) needs no weight and still
  applies.

```julia
using SpecialFunctions: loggamma
# the built-in gamma kernel, as a template (Chen's modified f₂ shape):
mygamma = AsymmetricKernel(
    (c, x, b) -> c > 0 ? (ρ = x < 2b ? (x / b)^2 / 4 + 1 : x / b;
                          exp((ρ - 1) * log(c) - c / b - ρ * log(b) - loggamma(ρ))) : 0.0,
    (0.0, Inf),
    x -> x)                       # the AMISE weight w(x) = x, for rule = :plugin
density(summary, 1.0; kernel = mygamma)
```

References:

- S. X. Chen (1999), [*Beta kernel estimators for density functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
- S. X. Chen (2000), [*Probability density function estimation using gamma kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3).
"""
struct AsymmetricKernel{F,W}
  pdf::F                            # (c, x, b) → a density in the data argument c
  support::Tuple{Float64,Float64}   # the data domain the kernel covers
  weight::W                         # the AMISE local-scale weight w(x) (x for the gamma kernel, x(1−x)
                                    # for the beta): both the leading bias ½b·w(x)f''(x) and the variance
                                    # ∝ w(x)^{-1/2} carry it. Used by the §4.5 plug-in bandwidth.
                                    # `nothing` (the 2-arg constructor) ⇒ rule=:plugin is unavailable.
end

# A custom kernel that only needs the rule-of-thumb bandwidth can omit the plug-in weight.
AsymmetricKernel(pdf, support) = AsymmetricKernel(pdf, support, nothing)

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

The **gamma kernel** (Chen 2000) for a density on the half-line ``[0,\\infty)`` — a varying-shape
*asymmetric* kernel: at the query point ``x`` it is the modified (bias-reduced) ``\\mathrm{Gamma}(\\rho(x),\\,b)``
density in the data argument (Chen's ``f_2``: shape ``\\rho(x) = x/b`` in the interior — so the kernel mean
is exactly ``x``, removing the ``f'`` bias of the naive ``x/b+1`` form — and ``(x/b)^2/4 + 1`` within ``2b``
of the edge), so its support matches the data and **no mass leaks past 0**. Passed as
`kernel = gamma_kernel`, [`density`](@ref GBGStreamingStatistics.DensityEstimation.density) returns the
boundary-bias-free estimate ``\\hat f(x) = \\tfrac1n \\sum_j n_j\\,\\mathrm{Gamma}(c_j;\\,\\rho(x),\\,b)``
(``O(b)`` bias uniformly, including at the boundary where a symmetric kernel is inconsistent; the variance
self-adapts, tightening toward the edge). Use it for non-negative quantities such as droplet radii — an
alternative to (and usually sharper than) the `bounds` renormalization. The bandwidth defaults to
``\\hat\\sigma\\,n^{-2/5}`` (the gamma-kernel rate; pass `h` to override, or the data-driven Chen plug-in
[`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth)`(summary; rule = :plugin, kernel = gamma_kernel)`,
``\\approx 1.4\\text{–}1.7\\times`` lower ISE on skewed data). `degree`, `adaptive`, and `bounds` do not
apply — the kernel is its own degree-0, boundary-free estimator.

Reference: S. X. Chen (2000), [*Probability density function estimation using gamma
kernels*](https://doi.org/10.1023/A:1004165218295), Ann. Inst. Statist. Math. 52(3).
"""
const gamma_kernel = AsymmetricKernel(_gamma_pdf, (0.0, Inf), x -> x)


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

The **beta kernel** (Chen 1999) for a density on the unit interval ``[0,1]`` — the bounded-interval
analog of [`gamma_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gamma_kernel): at ``x`` it is
the modified (bias-reduced) ``\\mathrm{Beta}`` density in the data argument (Chen's ``f_2``: interior
shapes ``(x/b,\\,(1-x)/b)``, with a ``\\sqrt{\\,\\cdot\\,}`` correction within ``2b`` of an edge),
boundary-bias-free at **both** edges. Pass `kernel = beta_kernel` with data scaled to ``[0,1]``; the bandwidth default, the knob
restrictions, and the ``n^{-2/5}`` rate are as for `gamma_kernel`.

Reference: S. X. Chen (1999), [*Beta kernel estimators for density
functions*](https://doi.org/10.1016/S0167-9473(99)00010-9), Comput. Statist. Data Anal. 31(2).
"""
const beta_kernel = AsymmetricKernel(_beta_pdf, (0.0, 1.0), x -> x * (1 - x))


# --- parametric-start (semiparametric) density (Hjort & Glad 1995) ---
# f̂(x) = p(x; θ̂)·r̂(x): fit a parametric "start" p(·; θ̂) — its parameters by method of moments from the
# recovered sample mean and variance, so θ̂ is itself a post-hoc functional of the summary — then KDE the
# *correction* r̂(x) = (1/n) Σᵢ K_h(x−xᵢ)/p(xᵢ; θ̂), the start's local relative error. The estimate is
# exact wherever the start is correct (the correction is then ≡ 1, so f̂ ≈ the parametric fit) and degrades
# only to the bias of a KDE-of-the-ratio where it is not — a tail-friendly bias reduction. Binned form:
# r̂(x) = (1/nh) Σⱼ (nⱼ/p(cⱼ; θ̂)) K((x−cⱼ)/h), a reduction over the bin counts and mergeable (θ̂ merges
# through the moments; the reweighted counts add). A start carries its *log*-density (stable: the
# correction weight is exp(−logpdf)) and its support; off-support points contribute nothing.
"""
    ParametricStart(fit, logpdf, support)

A **parametric start** for the Hjort–Glad semiparametric density (Hjort & Glad 1995). Passed as
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density)`(…; start = …)`, it multiplies the
moment-fitted model ``p(\\cdot;\\hat\\theta)`` by a kernel estimate of the **correction** ``f_X/p`` —
near-exact where the model fits (the bias is driven by ``(f_X/p)''``, not ``f_X''``), degrading only to a
kernel-of-the-ratio bias where it does not. The built-ins are
[`normal_start`](@ref GBGStreamingStatistics.DensityEstimation.normal_start),
[`lognormal_start`](@ref GBGStreamingStatistics.DensityEstimation.lognormal_start),
[`gamma_start`](@ref GBGStreamingStatistics.DensityEstimation.gamma_start),
[`exponential_start`](@ref GBGStreamingStatistics.DensityEstimation.exponential_start), and
[`inverse_gaussian_start`](@ref GBGStreamingStatistics.DensityEstimation.inverse_gaussian_start);
construct one directly to plug in a custom parametric model. The three fields are the construction contract:

- `fit(mean, variance)` — fit the model parameters from the **recovered moments** by method of moments,
  returning a parameter tuple. Fitting from the mean and variance (not the raw data) is what keeps the whole
  estimate **mergeable**: ``\\hat\\theta`` is itself a post-hoc functional of the summary.
- `logpdf(x, params)` — the model log-density ``\\log p(x;\\hat\\theta)``, returning ``-\\infty`` off the
  support. `density` reweights each bin by ``1/p(c_j;\\hat\\theta)`` through this log form (stable in the
  tails, where ``p`` is tiny).
- `support` — the domain the model is positive on, as a `(lower, upper)` tuple (e.g. `(0.0, Inf)`).

```jldoctest
julia> rayleigh = ParametricStart(
           (mean, variance) -> (mean * sqrt(2 / π),),                            # σ from the mean
           (x, p) -> x > 0 ? log(x) - 2 * log(p[1]) - x^2 / (2 * p[1]^2) : -Inf,  # Rayleigh log-density
           (0.0, Inf));

julia> rayleigh isa ParametricStart
true
```

Reference: N. L. Hjort & I. K. Glad (1995), [*Nonparametric density estimation with a parametric
start*](https://doi.org/10.1214/aos/1176324627), Ann. Statist. 23(3).
"""
struct ParametricStart{F,G}
  fit::F                            # (mean, variance) → params (a tuple), by method of moments
  logpdf::G                         # (x, params) → log p(x; params); −Inf off the start's support
  support::Tuple{Float64,Float64}   # where the start is positive (the estimand's domain)
end


"""
    normal_start

A Gaussian (normal) [parametric start](@ref GBGStreamingStatistics.DensityEstimation.density) on
``(-\\infty,\\infty)``: ``p(x;\\hat\\mu,\\hat\\sigma^2)`` with ``\\hat\\mu``, ``\\hat\\sigma^2`` the
recovered sample mean and variance. Pass `start = normal_start` to
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density); the estimate is then near-exact on
genuinely Gaussian data and reduces to a bias-corrected KDE otherwise.
"""
const normal_start = ParametricStart(
  (mean, variance) -> (mean, variance),
  (x, p) -> -0.5 * (log(2π) + log(p[2]) + (x - p[1])^2 / p[2]),
  (-Inf, Inf))


"""
    lognormal_start

A log-normal [parametric start](@ref GBGStreamingStatistics.DensityEstimation.density) on
``(0,\\infty)``, parameters by method of moments (``s^2 = \\log(1+\\hat\\sigma^2/\\hat\\mu^2)``,
``m = \\log\\hat\\mu - s^2/2``). A natural start for a **positively skewed** quantity such as a droplet
radius; pass `start = lognormal_start`.
"""
const lognormal_start = ParametricStart(
  (mean, variance) -> (log(mean) - log1p(variance / mean^2) / 2, log1p(variance / mean^2)),
  (x, p) -> x > 0 ? -log(x) - 0.5 * (log(2π) + log(p[2]) + (log(x) - p[1])^2 / p[2]) : -Inf,
  (0.0, Inf))


"""
    gamma_start

A gamma [parametric start](@ref GBGStreamingStatistics.DensityEstimation.density) on ``(0,\\infty)``,
shape ``k = \\hat\\mu^2/\\hat\\sigma^2`` and scale ``\\theta = \\hat\\sigma^2/\\hat\\mu`` by method of
moments. For a non-negative quantity; pass `start = gamma_start`.
"""
const gamma_start = ParametricStart(
  (mean, variance) -> (mean^2 / variance, variance / mean),
  (x, p) -> x > 0 ? (p[1] - 1) * log(x) - x / p[2] - p[1] * log(p[2]) - loggamma(p[1]) : -Inf,
  (0.0, Inf))


"""
    exponential_start

An exponential [parametric start](@ref GBGStreamingStatistics.DensityEstimation.density) on
``[0,\\infty)``, rate ``\\hat\\lambda = 1/\\hat\\mu`` (a one-parameter start fit from the mean alone).
Pass `start = exponential_start`; pairs with the `degree = 1` log-linear fit as the two
exponential-tail tools.
"""
const exponential_start = ParametricStart(
  (mean, variance) -> (1 / mean,),
  (x, p) -> x ≥ 0 ? log(p[1]) - p[1] * x : -Inf,
  (0.0, Inf))


"""
    inverse_gaussian_start

An inverse-Gaussian [parametric start](@ref GBGStreamingStatistics.DensityEstimation.density) on
``(0,\\infty)``, mean ``\\hat\\mu`` and shape ``\\hat\\lambda = \\hat\\mu^3/\\hat\\sigma^2`` by method of
moments — a positive start with an exponentially decaying right tail (lighter than the log-normal's).
Pass `start = inverse_gaussian_start`.
"""
const inverse_gaussian_start = ParametricStart(
  (mean, variance) -> (mean, mean^3 / variance),
  (x, p) -> x > 0 ? 0.5 * (log(p[2]) - log(2π) - 3 * log(x)) - p[2] * (x - p[1])^2 / (2 * p[1]^2 * x) : -Inf,
  (0.0, Inf))


# Validate the recovered moments against the start's support and fit θ̂ (method of moments).
function _start_params(summary::Histogram1D, start::ParametricStart)
  mean = _sample_mean(summary)
  variance = _sample_variance(summary)
  variance > 0 ||
    throw(ArgumentError("a parametric `start` needs a positive sample variance (degenerate summary; got variance = $variance)"))
  (isfinite(start.support[1]) && mean ≤ start.support[1]) &&
    throw(ArgumentError("a positive-support parametric `start` needs the sample mean inside its support (got mean = $mean, support lower bound = $(start.support[1]))"))
  return start.fit(mean, variance)
end


# Reject the knobs that do not combine with a parametric start (it is the degree-0, fixed-h, symmetric-
# kernel local-constant estimator of Hjort & Glad; `normalize` is its renormalization, `bounds` its rival).
function _check_start(kernel, degree, adaptive, bounds)
  kernel isa AsymmetricKernel &&
    throw(ArgumentError("a parametric `start` pairs with a symmetric kernel; the asymmetric `gamma_kernel`/`beta_kernel` are already boundary-free — drop `start`"))
  (degree == 0 && !adaptive && bounds === nothing) ||
    throw(ArgumentError("a parametric `start` is the degree-0, fixed-bandwidth estimator without `bounds`; `degree ≥ 1`, `adaptive`, and `bounds` do not combine with it (use `normalize` to enforce ∫f̂ = 1)"))
end


# Hjort–Glad estimate at x: f̂(x) = p(x; θ̂)·(1/nh) Σⱼ (nⱼ/p(cⱼ; θ̂)) K((x−cⱼ)/h). Off the start's support
# (p = 0) the estimate is 0, and bins off the support carry no correction weight.
function _parametric_start_eval(counts, centers, n, x, h, kernel, start, params)
  logpx = start.logpdf(x, params)
  isfinite(logpx) || return 0.0
  invh = inv(h)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    logpj = start.logpdf(centers[j], params)
    isfinite(logpj) || continue
    acc += cj * exp(-logpj) * kernel((x - centers[j]) * invh)
  end
  return exp(logpx) * acc / (n * h)
end


# --- exact renormalization: divide the estimate by ∫f̂ over the summary's range ∩ the estimator's support
# so ∫f̂ = 1 exactly. The symmetric KDE already integrates to ≈1; the asymmetric kernels and a parametric
# start do not, so `normalize` is their finishing step. A post-hoc functional (off the hot path). ---
# The estimand's support: the asymmetric kernel's domain, the `bounds`, and/or the start's support, met.
function _density_support(kernel, bounds, start)
  lower, upper = -Inf, Inf
  kernel isa AsymmetricKernel && ((lower, upper) = kernel.support)
  bounds === nothing || ((lower, upper) = (max(lower, bounds[1]), min(upper, bounds[2])))
  start === nothing || ((lower, upper) = (max(lower, start.support[1]), min(upper, start.support[2])))
  return lower, upper
end


# Z = ∫ f̂ over [summary.lower, summary.upper] ∩ support, by the trapezoidal rule on a fine grid.
function _normalization_constant(summary::Histogram1D; degree, h, adaptive, kernel, bounds, start, binning)
  slower, supper = _density_support(kernel, bounds, start)
  lower = max(float(summary.lower), slower)
  upper = min(float(summary.upper), supper)
  npts = clamp(4 * _nbins(summary), 512, 8192)
  grid = collect(range(lower, upper; length = npts))
  vals = density(summary, grid; degree = degree, h = h, adaptive = adaptive,
                 kernel = kernel, bounds = bounds, start = start, normalize = false, binning = binning)
  Δ = (upper - lower) / (npts - 1)
  return (sum(vals) - 0.5 * (vals[1] + vals[end])) * Δ
end


# --- default bandwidth: the kernel- and degree-aware normal-reference rule (§4.5) ---
# Kernel constants come from numerical quadrature of K, so any user-supplied `kernel` works: integrate
# g over [−L, L] by the trapezoidal rule on a fixed fine grid. The shipped kernels sit well inside it
# (Gaussian ≈ ±8, Epanechnikov ±1); this runs once per `density`/`bandwidth` call, off the hot path.
function _quad(g, L::Float64 = 16.0, npts::Int = 32_769)
  Δ = 2L / (npts - 1)
  acc = 0.5 * (g(-L) + g(L))
  @inbounds for i in 1:npts-2
    acc += g(-L + i * Δ)
  end
  return acc * Δ
end


_doublefactorial(k::Integer) = k ≤ 1 ? 1.0 : k * _doublefactorial(k - 2)


# The canonical kernel×degree factor Wikk(K, deg) and the leading bias order p — the same grouping as
# locfit's `Wikk` (Loader 1996). Degree 0,1 share the second-order kernel: Wikk = R(K)/μ₂(K)², p = 2.
# Degree 2,3 use the fourth-order *equivalent* kernel K⋆ = K·(n₄ − v²n₂)/(n₀n₄ − n₂²): Wikk =
# (4!²/8)·R(K⋆)/μ₄(K⋆)², p = 4 — with kernel moments nₖ = ∫vᵏK. Higher degree is a recorded coverage gap.
function _wikk(kernel, degree::Integer)
  n0 = _quad(kernel)
  n2 = _quad(v -> v^2 * kernel(v))
  if degree == 0 || degree == 1
    return _quad(v -> kernel(v)^2) / n2^2, 2
  elseif degree == 2 || degree == 3
    n4 = _quad(v -> v^4 * kernel(v))
    determinant = n0 * n4 - n2^2
    Kstar = v -> kernel(v) * (n4 - v^2 * n2) / determinant
    return (factorial(4)^2 / 8) * _quad(v -> Kstar(v)^2) / _quad(v -> v^4 * Kstar(v))^2, 4
  end
  throw(ArgumentError("normal-reference bandwidth implemented for degree ≤ 3 (got degree = $degree)"))
end


# The kernel- and degree-aware normal-reference (Gaussian-reference) bandwidth — the default `h`:
#   h = C(K, deg)·σ̂·n^(−1/(2p+1)),   C = [Wikk(K, deg)·2^(p+1)·√π / (2p−1)!!]^(1/(2p+1)),
# recovered from σ̂ and n. Reduces to Silverman (4/3)^(1/5)·σ̂·n^(−1/5) for the Gaussian kernel at
# degree 0/1; degree 2 is a higher-order (O(h⁴)) fit, so it widens to the n^(−1/9) rate. See `bandwidth`.
function _normal_reference_h(summary::Histogram1D, kernel, degree::Integer)
  wikk, p = _wikk(kernel, degree)
  C = (wikk * 2.0^(p + 1) * sqrt(π) / _doublefactorial(2p - 1))^(1 / (2p + 1))
  return C * _sample_std(summary) * float(_total(summary))^(-1 / (2p + 1))
end


# Fixed-bandwidth binned KDE at x (single h) — the degree-0 fixed estimator and the adaptive pilot.
function _kde_h(counts, centers, n, x, h, kernel)
  invh = inv(h)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel((x - centers[j]) * invh)
  end
  return acc / (n * h)
end


# Abramson: per-bin (sample-point) bandwidths ``h_j = h (\tilde f(c_j)/g)^{-1/2}`` and the geometric
# mean ``g`` of the fixed-h pilot ``\tilde f`` over the sample (occupied bins, weighted by counts). The
# ``h_j`` drive the degree-0 adaptive KDE; ``g`` drives the degree-≥1 balloon bandwidth. O(nbins²).
function _abramson(counts, centers, n, h, kernel)
  nb = length(centers)
  pilot = zeros(Float64, nb)
  @inbounds for j in 1:nb
    iszero(counts[j]) && continue
    pilot[j] = _kde_h(counts, centers, n, centers[j], h, kernel)
  end
  logg = 0.0
  @inbounds for j in 1:nb
    (counts[j] > 0 && pilot[j] > 0) || continue
    logg += counts[j] * log(pilot[j])
  end
  g = exp(logg / n)
  hj = [pilot[j] > 0 ? h * sqrt(g / pilot[j]) : h for j in 1:nb]
  return hj, g
end


# Degree-0 adaptive KDE at x with per-bin (sample-point) bandwidths ``h_j``.
function _kde_adaptive(counts, centers, hj, n, x, kernel)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel((x - centers[j]) / hj[j]) / hj[j]
  end
  return acc / n
end


# Degree ≥ 1 — local-log-polynomial density by binned **local Poisson likelihood** (Loader 1996), at a
# single point-wise (balloon) bandwidth ``h_x``. Model the bin counts ``n_j`` as Poisson with
# ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l=0}^{d}\beta_l (c_j-x)^l`` and weights ``K((c_j-x)/h_x)``;
# Newton/IRLS (concave) gives ``\hat f(x) = e^{\beta_0}``. Warm-started from the KDE value `f0`; falls
# back to it where the local design can't support degree d (singular info matrix) or the fit diverges.
function _local_log(counts, centers, hx, n, b, x, degree, kernel, f0)
  (f0 > 0 && isfinite(f0)) || return f0
  invhx = inv(hx)
  p = degree + 1
  β = zeros(p)
  β[1] = log(f0)
  nb_off = n * b
  grad = zeros(p)
  H = zeros(p, p)
  pw = zeros(p)
  for _ in 1:25
    fill!(grad, 0.0)
    fill!(H, 0.0)
    @inbounds for j in eachindex(centers)
      d = centers[j] - x
      w = kernel(d * invhx)
      w < 1e-10 && continue
      pw[1] = 1.0
      η = β[1]
      acc = 1.0
      for l in 2:p
        acc *= d
        pw[l] = acc
        η += β[l] * acc
      end
      μ = nb_off * exp(η)
      r = counts[j] - μ
      for l in 1:p
        grad[l] += w * r * pw[l]
        wμpl = w * μ * pw[l]
        for m in 1:p
          H[l, m] += wμpl * pw[m]
        end
      end
    end
    δ = try
      H \ grad
    catch
      return f0
    end
    all(isfinite, δ) || return f0
    β .+= δ
    maximum(abs, δ) < 1e-8 && break
  end
  fx = exp(β[1])
  # a bounded refinement of the KDE: a fit that departs by more than ×5 is an extrapolation artifact
  # (a log-polynomial run wild in a sparse / wide-window tail), so fall back to the KDE there.
  return (isfinite(fx) && f0 / 5 ≤ fx ≤ 5 * f0) ? fx : f0
end


# The point-wise (balloon) bandwidth at x for the local-log fit: fixed `h`, or the Abramson factor
# ``h (\tilde f(x)/g)^{-1/2}`` evaluated *at the query point* (so the fit window stays local). The factor
# is clipped at 5 so the window does not blow up in the deep tail (where ``\tilde f(x)\to 0``) and the
# local polynomial extrapolate upward.
@inline _balloon(h, g, f0, adaptive) = (adaptive && f0 > 0) ? h * min(sqrt(g / f0), 5.0) : float(h)


# --- boundary renormalization (Jones 1993): divide the degree-0 KDE by the kernel mass that lands inside
# the support, so a density on a half-line / interval is not biased down at the edge (e.g. droplet radii
# on [0, ∞), whose KDE otherwise drops to ≈ half the true value at the boundary). ---
# The kernel's CDF (mass to the left of t), by quadrature of K — reuses the fixed grid of `_quad`.
_kernel_cdf(kernel, t) = _quad(s -> s ≤ t ? kernel(s) : 0.0)


# The fraction of a kernel centred at x (bandwidth h) inside the support [a, b]: M(x) = ∫_a^b K_h(x−u)du =
# CDF_K((x−a)/h) − CDF_K((x−b)/h)  (a = −∞ ⇒ the first term is 1; b = +∞ ⇒ the second is 0).
function _kernel_mass(kernel, x, h, a, b)
  upper = isfinite(a) ? _kernel_cdf(kernel, (x - a) / h) : 1.0
  lower = isfinite(b) ? _kernel_cdf(kernel, (x - b) / h) : 0.0
  return upper - lower
end


# Asymmetric-kernel estimate f̂(x) = (1/n) Σ_j n_j K(c_j; x, b) — a reduction over the bin counts.
function _asymmetric_eval(counts, centers, n, kernel, x, b)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel(centers[j], x, b)
  end
  return acc / n
end

# Default bandwidth for an asymmetric kernel: the gamma/beta-kernel rate b ∝ n^(−2/5) (Chen 1999/2000),
# scaled by σ̂ so the kernel's width at a typical point matches the symmetric optimum.
_asymmetric_default_b(summary::Histogram1D) = _sample_std(summary) * float(_total(summary))^(-2 / 5)


# Each bin's representative location for the post-hoc estimators, set by `binning`. `:simple` — the
# geometric bin center c_j (the histogram's mass-at-center approximation, an O(Δ) discretization of the
# exact KDE). `:linear` — the recovered within-bin centroid ĉ_j = c_j + S_1[j]/n_j from the stored first
# power-sum, which matches the bin's own first moment, so the within-bin displacement averages to zero and
# the binned KDE's discretization error drops to O(Δ²) (the centroid analog of linear binning). Every eval
# helper is agnostic to which — it just reads `centers`. Empty bins keep the geometric center (no mass).
function _bin_locations(summary::Histogram1D, counts, binning::Symbol)
  nb = length(counts)
  if binning === :simple
    return [_bincenter(summary, j) for j in 1:nb]
  elseif binning === :linear
    S1 = Array(summary.powersums[1, :])      # the first power-sum S_1 about each bin center
    return [iszero(counts[j]) ? _bincenter(summary, j) : _bincenter(summary, j) + S1[j] / counts[j]
            for j in 1:nb]
  end
  throw(ArgumentError("binning must be :simple or :linear (got :$(binning))"))
end


"""
    density(summary::Histogram1D, x; degree = 0, h = nothing, adaptive = false, kernel = gaussian_kernel, bounds = nothing, start = nothing, normalize = false)   # ≡ f_X

**Estimand.** the marginal (Lebesgue) density ``f_X = \\mathrm{d}\\mu_X/\\mathrm{d}\\lambda`` of the
conditioning variable ``X``, defined by ``P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda``.

**Estimator.** a binned local smoother of the counts ``n_j`` at centers ``c_j``, with these knobs:

- **`degree`** — the local-polynomial order. `degree = 0` (default) is the binned KDE / WARPing
  ``\\hat f(x) = \\frac1n \\sum_j (n_j/h_j)\\, K((x-c_j)/h_j)`` (``\\int K = 1``). `degree` ``\\ge 1`` is the
  **local-log-polynomial** (local-likelihood) estimator (Loader 1996): fit ``\\log f`` by a
  degree-`degree` polynomial through a binned **local Poisson** likelihood and report
  ``\\hat f(x) = e^{\\hat\\beta_0}``. `degree = 1` (log-linear) is **exact for an exponential density**,
  so it tracks the exponential / stretched-exponential **tails** a flat KDE smears; `degree = 2`
  (log-quadratic) is exact for a Gaussian. It falls back to the KDE where the local fit is degenerate or
  runs away in a sparse tail.
- **`adaptive`** — the bandwidth. `false` (default) uses the fixed `h`; `true` uses the **Abramson**
  variable bandwidth — the kernel width scales by ``(\\tilde f/g)^{-1/2}`` (a pilot ``\\tilde f`` at `h`,
  ``g`` its geometric mean), widening in the sparse tails: far less tail variance (Abramson 1982). `h`
  is the fixed bandwidth or the pilot; it **defaults to the kernel- and degree-aware normal-reference
  rule** ([`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth), §4.5 — Silverman's
  ``(4/3)^{1/5}\\hat\\sigma\\,n^{-1/5}`` for the Gaussian kernel at degree 0/1, wider at the
  ``n^{-1/9}`` rate for degree 2). (The KDE uses the sample-point ``h_j`` per bin; the local-log fit
  uses one balloon ``h(x)`` at the query.)
- **`bounds`** — the support of the density, passed as the tuple `(lower, upper)` (e.g. `(0.0, Inf)` for
  a positive quantity such as a droplet radius), or `nothing` (default, unbounded). When given, the
  degree-0 KDE is **renormalized** by the kernel mass inside the support: ``\\hat f(x) \\big/ \\int_a^b K_h(x-u)\\,\\mathrm{d}u`` (Jones 1993). This corrects the boundary bias that otherwise drops the estimate
  to ``\\approx`` half the true value at the edge. (`degree` ``\\ge 1`` is already boundary-adaptive, so
  `bounds` is a degree-0 option.)
- **`start`** — a **parametric start** (Hjort & Glad 1995): the semiparametric estimator
  ``\\hat f(x) = p(x;\\hat\\theta)\\,\\frac1{nh}\\sum_j \\frac{n_j}{p(c_j;\\hat\\theta)}\\,K((x-c_j)/h)``,
  which multiplies a fitted parametric model ``p(\\cdot;\\hat\\theta)`` by a KDE of its *correction*
  (the local relative error). ``\\hat\\theta`` is fit by method of moments from the recovered
  ``\\hat\\mu``, ``\\hat\\sigma^2``, so the start is itself a functional of the summary. The estimate is
  **near-exact wherever the start is correct** (the correction is then ``\\equiv 1``) and degrades only to
  a bias-corrected KDE where it is not — a tail-friendly bias reduction. Pass one of
  [`normal_start`](@ref GBGStreamingStatistics.DensityEstimation.normal_start),
  [`lognormal_start`](@ref GBGStreamingStatistics.DensityEstimation.lognormal_start),
  [`gamma_start`](@ref GBGStreamingStatistics.DensityEstimation.gamma_start),
  [`exponential_start`](@ref GBGStreamingStatistics.DensityEstimation.exponential_start),
  [`inverse_gaussian_start`](@ref GBGStreamingStatistics.DensityEstimation.inverse_gaussian_start). It is
  a degree-0, fixed-`h`, symmetric-kernel estimator (`degree` ``\\ge 1``, `adaptive`, `bounds`, and the
  asymmetric kernels do not combine with it; use `normalize` for ``\\int\\hat f = 1``).
- **`normalize`** — when `true`, divide the estimate by ``\\int\\hat f`` over the summary's range
  (intersected with the estimator's support), so ``\\int\\hat f = 1`` **exactly**. The symmetric KDE
  already integrates to ``\\approx 1``; the asymmetric kernels and a parametric start do not, so this is
  their finishing step.
- **`binning`** — where each bin's mass sits for the post-hoc evaluation. `:linear` (the default) places
  it at the recovered within-bin **centroid** ``c_j + S_1^{(j)}/n_j`` (from the stored first power-sum
  ``S_1``), which matches the bin's own first moment, so the binned estimate converges to the exact KDE at
  ``O(\\Delta^2)`` in the bin width ``\\Delta``; `:simple` places it at the geometric bin center ``c_j``,
  the ``O(\\Delta)`` mass-at-center approximation. The centroid is read from data already in the summary, so
  the higher-order accuracy is free; it applies to every path (KDE, local-log, start, asymmetric kernel).

The four corners (fixed/variable ``\\times`` KDE/local-log) and everything between are one call. `kernel`
is any ``K`` with ``\\int K = 1`` ([`gaussian_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gaussian_kernel)
default, [`epanechnikov_kernel`](@ref GBGStreamingStatistics.DensityEstimation.epanechnikov_kernel)); for
**bounded / positive support** pass the asymmetric
[`gamma_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gamma_kernel) (half-line) or
[`beta_kernel`](@ref GBGStreamingStatistics.DensityEstimation.beta_kernel) (unit interval), which are
boundary-bias-free by construction. A vector of query points returns a vector. The FFT amortization of the `degree = 0`, fixed-`h` convolution
over a bandwidth grid is forthcoming (`CLAUDE.md` roadmap).

```jldoctest
julia> h = Histogram1D(-1.0, 1.0, 2);

julia> observe!(h, [0.4]);          # one sample at 0.4

julia> round(density(h, 0.4; h = 0.2); digits = 4)   # linear binning centers the bump at 0.4 (the centroid)
1.9947
```

References:

- I. S. Abramson (1982), [*On bandwidth variation in kernel estimates—a square root law*](https://doi.org/10.1214/aos/1176345986), Ann. Statist. 10(4) (variable bandwidth).
- C. R. Loader (1996), [*Local likelihood density estimation*](https://doi.org/10.1214/aos/1032298287), Ann. Statist. 24(4) (the local-log-polynomial estimator).
- M. C. Jones (1993), [*Simple boundary correction for kernel density estimation*](https://doi.org/10.1007/BF00147776), Statist. Comput. 3(3) (the `bounds` renormalization).
- N. L. Hjort & I. K. Glad (1995), [*Nonparametric density estimation with a parametric start*](https://doi.org/10.1214/aos/1176324627), Ann. Statist. 23(3) (the `start` estimator).
"""
function density(summary::Histogram1D, x::Real; degree::Integer = 0, h = nothing,
                 adaptive::Bool = false, kernel = gaussian_kernel, bounds = nothing,
                 start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  if normalize
    Z = _normalization_constant(summary; degree = degree, h = h, adaptive = adaptive,
                                kernel = kernel, bounds = bounds, start = start, binning = binning)
    fx = density(summary, x; degree = degree, h = h, adaptive = adaptive, kernel = kernel,
                 bounds = bounds, start = start, normalize = false, binning = binning)
    return Z > 0 ? fx / Z : fx
  end
  if start !== nothing
    _check_start(kernel, degree, adaptive, bounds)
    n = _total(summary)
    iszero(n) && return 0.0
    hh = h === nothing ? _normal_reference_h(summary, kernel, 0) : float(h)
    hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    params = _start_params(summary, start)
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return _parametric_start_eval(counts, centers, n, x, hh, kernel, start, params)
  end
  if kernel isa AsymmetricKernel
    (degree == 0 && !adaptive && bounds === nothing) ||
      throw(ArgumentError("the asymmetric kernel is its own degree-0, boundary-free estimator; " *
                          "`degree`, `adaptive`, and `bounds` do not apply"))
    n = _total(summary)
    iszero(n) && return 0.0
    bb = h === nothing ? _asymmetric_default_b(summary) : float(h)
    bb > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return _asymmetric_eval(counts, centers, n, kernel, x, bb)
  end
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  (bounds === nothing || degree == 0) ||
    throw(ArgumentError("`bounds` (boundary renormalization) applies to the degree-0 KDE; " *
                        "degree ≥ 1 (local-log) is already boundary-adaptive (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return 0.0
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = _bin_locations(summary, counts, binning)
  if degree == 0 && !adaptive
    fx = _kde_h(counts, centers, n, x, hh, kernel)
    return bounds === nothing ? fx : fx / _kernel_mass(kernel, x, hh, bounds[1], bounds[2])
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    fx = _kde_adaptive(counts, centers, hj, n, x, kernel)
    return bounds === nothing ? fx : fx / _kernel_mass(kernel, x, hh, bounds[1], bounds[2])
  end
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  f0 = _kde_h(counts, centers, n, x, hh, kernel)
  return _local_log(counts, centers, _balloon(hh, g, f0, adaptive), n, inv(summary.scale),
                    x, degree, kernel, f0)
end

function density(summary::Histogram1D, xs::AbstractVector; degree::Integer = 0, h = nothing,
                 adaptive::Bool = false, kernel = gaussian_kernel, bounds = nothing,
                 start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  if normalize
    Z = _normalization_constant(summary; degree = degree, h = h, adaptive = adaptive,
                                kernel = kernel, bounds = bounds, start = start, binning = binning)
    fxs = density(summary, xs; degree = degree, h = h, adaptive = adaptive, kernel = kernel,
                  bounds = bounds, start = start, normalize = false, binning = binning)
    return Z > 0 ? fxs ./ Z : fxs
  end
  if start !== nothing
    _check_start(kernel, degree, adaptive, bounds)
    n = _total(summary)
    iszero(n) && return zeros(Float64, length(xs))
    hh = h === nothing ? _normal_reference_h(summary, kernel, 0) : float(h)
    hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    params = _start_params(summary, start)
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return [_parametric_start_eval(counts, centers, n, x, hh, kernel, start, params) for x in xs]
  end
  if kernel isa AsymmetricKernel
    (degree == 0 && !adaptive && bounds === nothing) ||
      throw(ArgumentError("the asymmetric kernel is its own degree-0, boundary-free estimator; " *
                          "`degree`, `adaptive`, and `bounds` do not apply"))
    n = _total(summary)
    iszero(n) && return zeros(Float64, length(xs))
    bb = h === nothing ? _asymmetric_default_b(summary) : float(h)
    bb > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return [_asymmetric_eval(counts, centers, n, kernel, x, bb) for x in xs]
  end
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  (bounds === nothing || degree == 0) ||
    throw(ArgumentError("`bounds` (boundary renormalization) applies to the degree-0 KDE; " *
                        "degree ≥ 1 (local-log) is already boundary-adaptive (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = _bin_locations(summary, counts, binning)
  if degree == 0 && !adaptive
    fxs = [_kde_h(counts, centers, n, x, hh, kernel) for x in xs]
    return bounds === nothing ? fxs : [f / _kernel_mass(kernel, x, hh, bounds[1], bounds[2]) for (f, x) in zip(fxs, xs)]
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    fxs = [_kde_adaptive(counts, centers, hj, n, x, kernel) for x in xs]
    return bounds === nothing ? fxs : [f / _kernel_mass(kernel, x, hh, bounds[1], bounds[2]) for (f, x) in zip(fxs, xs)]
  end
  b = inv(summary.scale)
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  return [begin
            f0 = _kde_h(counts, centers, n, x, hh, kernel)
            _local_log(counts, centers, _balloon(hh, g, f0, adaptive), n, b, x, degree, kernel, f0)
          end for x in xs]
end


"""
    density_standard_error(summary::Histogram1D, x; h = nothing, kernel = gaussian_kernel, degree = 0)   # ≡ SE[f̂(x)]

**Estimand.** the pointwise standard error of the density estimator
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density) at ``x`` — ``\\mathrm{SE}(x) =
\\sqrt{\\operatorname{Var}\\hat f_h(x)}`` with the kernel-density asymptotic variance

```math
\\operatorname{Var}\\hat f_h(x) \\approx \\frac{R(K)\\,f_X(x)}{n h}, \\qquad R(K) = \\int K^2.
```

**Estimator.** ``\\sqrt{R(K)\\,\\hat f_h(x)/(n h)}`` from the recovered ``\\hat f_h(x)``, the count
``n``, the bandwidth ``h`` (default the §4.5 normal reference), and ``R(K)`` by quadrature. The
pointwise **variability band** is ``\\hat f_h(x) \\pm z\\,\\mathrm{SE}(x)`` (``z = 1.96`` for 95\\%) — a
band about ``\\mathbb{E}\\hat f_h``, so it covers ``f_X`` where the ``O(h^2)`` smoothing bias is small
and narrows that coverage at a sharp mode or in a tail. `degree` 0/1 (the KDE variance; `degree` ``\\ge 2``
and the variable bandwidth are a coverage gap); a vector of query points returns a vector.

```jldoctest
julia> h = Histogram1D(-6.0, 6.0, 600);

julia> observe!(h, [-0.4, 0.1, 0.2, 0.9, 1.3]);

julia> density_standard_error(h, 0.5) > 0
true
```
"""
function density_standard_error(summary::Histogram1D, x::Real; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, binning::Symbol = :linear)
  (degree == 0 || degree == 1) ||
    throw(ArgumentError("density_standard_error is a degree-0/1 (KDE) variability band; degree ≥ 2 is a coverage gap (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return 0.0
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  fx = density(summary, x; h = hh, kernel = kernel, degree = degree, binning = binning)
  return sqrt(_quad(u -> kernel(u)^2) * max(fx, 0.0) / (n * hh))
end

function density_standard_error(summary::Histogram1D, xs::AbstractVector; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, binning::Symbol = :linear)
  (degree == 0 || degree == 1) ||
    throw(ArgumentError("density_standard_error is a degree-0/1 (KDE) variability band; degree ≥ 2 is a coverage gap (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  RK = _quad(u -> kernel(u)^2)
  fxs = density(summary, xs; h = hh, kernel = kernel, degree = degree, binning = binning)
  return [sqrt(RK * max(f, 0.0) / (n * hh)) for f in fxs]
end

end # module DensityEstimation
