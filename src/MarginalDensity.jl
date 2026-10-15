# MarginalDensity — §4.1: the marginal density f_X. A single `density` spans the local-log-polynomial
# degree (0 = KDE, ≥1 = local likelihood) crossed with a fixed or Abramson-variable bandwidth, plus
# pluggable kernels (the symmetric gaussian/epanechnikov and the bounded-support asymmetric gamma/beta),
# the Jones boundary `support` renormalization, the Hjort–Glad parametric `start`, exact `normalize`, and
# the `density_standard_error` band. Two stages: set the bandwidth, then fit. The bandwidth is
# *sample-point* (per-bin h_j) for the KDE and *balloon* (one h(x) at the query) for the local-log fit.
# Written estimand-first (`docs/ARCHITECTURE.md` §2, §4.1).

module MarginalDensity

using SpecialFunctions: loggamma

using ..Histograms: Histogram1D, _nbins, _bincenter, _total
using ..Expectations: _sample_mean, _sample_variance, _sample_std
using ..LocalPolynomial: _hankel_solve
using ..Kernels: gaussian_kernel, epanechnikov_kernel, _quad, _kernel_weight, AsymmetricKernel,
                 gamma_kernel, beta_kernel

export density, density_standard_error, density_with_standard_error,
       ParametricStart, normal_start, lognormal_start, gamma_start, exponential_start, inverse_gaussian_start


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
[`density(…; start = …)`](@ref GBGStreamingStatistics.MarginalDensity.density), it multiplies the
moment-fitted model ``p(\\cdot;\\hat\\theta)`` by a kernel estimate of the **correction** ``f_X/p`` —
near-exact where the model fits (the bias is driven by ``(f_X/p)''``, not ``f_X''``), degrading only to a
kernel-of-the-ratio bias where it does not. The built-ins are
[`normal_start`](@ref GBGStreamingStatistics.MarginalDensity.normal_start),
[`lognormal_start`](@ref GBGStreamingStatistics.MarginalDensity.lognormal_start),
[`gamma_start`](@ref GBGStreamingStatistics.MarginalDensity.gamma_start),
[`exponential_start`](@ref GBGStreamingStatistics.MarginalDensity.exponential_start), and
[`inverse_gaussian_start`](@ref GBGStreamingStatistics.MarginalDensity.inverse_gaussian_start);
construct one directly to plug in a custom parametric model.

# Fields

- `fit` — the moment fit, a callable `(mean, variance)` returning the parameter tuple.
- `logpdf` — the model log-density, a callable `(x, params)`.
- `support` — the domain the model is positive on, a `(lower, upper)` tuple.

# Arguments

The constructor arguments are the three fields, and they are the construction contract:

- `fit(mean, variance)` — fit the model parameters from the **recovered moments** by method of moments,
  returning a parameter tuple. Fitting from the mean and variance (not the raw data) is what keeps the whole
  estimate **mergeable**: ``\\hat\\theta`` is itself a post-hoc functional of the summary.
- `logpdf(x, params)` — the model log-density ``\\log p(x;\\hat\\theta)``, returning ``-\\infty`` off the
  support. `density` reweights each bin by ``1/p(\\ell_j;\\hat\\theta)`` at the bin location ``\\ell_j`` through
  this log form (stable in the tails, where ``p`` is tiny).
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

A Gaussian (normal) [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) on
``(-\\infty,\\infty)``: ``p(x;\\hat\\mu,\\hat\\sigma^2)`` with ``\\hat\\mu``, ``\\hat\\sigma^2`` the
recovered sample mean and variance. Pass `start = normal_start` to
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density).
"""
const normal_start = ParametricStart(
  (mean, variance) -> (mean, variance),
  (x, p) -> -0.5 * (log(2π) + log(p[2]) + (x - p[1])^2 / p[2]),
  (-Inf, Inf))


"""
    lognormal_start

A log-normal [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) on
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

A gamma [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) on ``(0,\\infty)``,
shape ``k = \\hat\\mu^2/\\hat\\sigma^2`` and scale ``\\theta = \\hat\\sigma^2/\\hat\\mu`` by method of
moments. For a non-negative quantity; pass `start = gamma_start`.
"""
const gamma_start = ParametricStart(
  (mean, variance) -> (mean^2 / variance, variance / mean),
  (x, p) -> x > 0 ? (p[1] - 1) * log(x) - x / p[2] - p[1] * log(p[2]) - loggamma(p[1]) : -Inf,
  (0.0, Inf))


"""
    exponential_start

An exponential [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) on
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

An inverse-Gaussian [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) on
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


# Reject the options that do not combine with a parametric start (it is the degree-0, fixed-h, symmetric-
# kernel local-constant estimator of Hjort & Glad; `normalize` is its renormalization, `support` its rival).
function _check_start(kernel, degree, adaptive, support)
  kernel isa AsymmetricKernel &&
    throw(ArgumentError("a parametric `start` pairs with a symmetric kernel; the asymmetric `gamma_kernel`/`beta_kernel` are already boundary-free — drop `start`"))
  (degree == 0 && !adaptive && support === nothing) ||
    throw(ArgumentError("a parametric `start` is the degree-0, fixed-bandwidth estimator without `support`; `degree ≥ 1`, `adaptive`, and `support` do not combine with it (use `normalize` to enforce ∫f̂ = 1)"))
end


# Hjort–Glad estimate at x: f̂(x) = p(x; θ̂)·(1/nh) Σⱼ (nⱼ/p(cⱼ; θ̂)) K((x−cⱼ)/h). Off the start's support
# (p = 0) the estimate is 0, and bins off the support carry no correction weight.
function _parametric_start_eval(counts, centers, n, x, h, kernel, start, params)
  logpx = start.logpdf(x, params)
  isfinite(logpx) || return 0.0
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    logpj = start.logpdf(centers[j], params)
    isfinite(logpj) || continue
    acc += cj * exp(-logpj) * _kernel_weight(centers[j], x, h, kernel)
  end
  return exp(logpx) * acc / n
end


# --- exact renormalization: divide the estimate by ∫f̂ over the summary's range ∩ the estimator's support
# so ∫f̂ = 1 exactly. The symmetric KDE already integrates to ≈1; the asymmetric kernels and a parametric
# start do not, so `normalize` is their finishing step. A post-hoc functional (off the hot path). ---
# The estimand's support: the asymmetric kernel's domain, the `support`, and/or the start's support, met.
function _density_support(kernel, support, start)
  lower, upper = -Inf, Inf
  kernel isa AsymmetricKernel && ((lower, upper) = kernel.support)
  support === nothing || ((lower, upper) = (max(lower, support[1]), min(upper, support[2])))
  start === nothing || ((lower, upper) = (max(lower, start.support[1]), min(upper, start.support[2])))
  return lower, upper
end


# Z = ∫ f̂ over [summary.lower, summary.upper] ∩ support, by the trapezoidal rule on a fine grid.
function _normalization_constant(summary::Histogram1D; degree, h, adaptive, kernel, support, start, binning)
  slower, supper = _density_support(kernel, support, start)
  lower = max(float(summary.lower), slower)
  upper = min(float(summary.upper), supper)
  npts = clamp(4 * _nbins(summary), 512, 8192)
  grid = collect(range(lower, upper; length = npts))
  vals = density(summary, grid; degree = degree, h = h, adaptive = adaptive,
                 kernel = kernel, support = support, start = start, normalize = false, binning = binning)
  Δ = (upper - lower) / (npts - 1)
  return (sum(vals) - 0.5 * (vals[1] + vals[end])) * Δ
end


# --- default bandwidth: the kernel- and degree-aware normal-reference rule (§4.5) ---
# Kernel constants (``R(K)``, ``\mu_2(K)``, ``\mu_4(K)``) come from the shared `Kernels._quad` quadrature
# of K, so any user-supplied `kernel` works with no hand-coded constants.


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


# Fixed-bandwidth binned degree-0 density at x (single h): f̂ = (1/n) Σ_j n_j w_j with w_j the shared
# `_kernel_weight` (the 1/h folded in) — so this one evaluator serves both the symmetric KDE and the
# asymmetric (gamma/beta) boundary kernel, and is the adaptive pilot. The degree-0 readout of the atom.
function _kde_h(counts, centers, n, x, h, kernel)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * _kernel_weight(centers[j], x, h, kernel)
  end
  return acc / n
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
function _kde_adaptive(counts, centers, n, x, hj, kernel)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * _kernel_weight(centers[j], x, hj[j], kernel)
  end
  return acc / n
end


# Degree ≥ 1 — local-log-polynomial density by binned **local Poisson likelihood** (Loader 1996), at a
# single point-wise (balloon) bandwidth ``h_x``. Model the bin counts ``n_j`` as Poisson with
# ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l=0}^{d}\beta_l (c_j-x)^l`` and weights ``K((c_j-x)/h_x)``;
# Newton/IRLS (concave) gives ``\hat f(x) = e^{\beta_0}``. Each Newton step is the shared
# kernel-weighted local-polynomial solve (`LocalPolynomial._hankel_solve` of the μ-weighted design
# moments) — a local GLM is the linear solve iterated (`docs/ARCHITECTURE.md` §11). Warm-started from
# the KDE value `f0`; falls back to it where the local design can't support degree d (singular) or the
# fit diverges.
function _local_log(counts, centers, n, x, hx, kernel, degree, b, f0)
  (f0 > 0 && isfinite(f0)) || return f0
  p = degree + 1
  β = zeros(p)
  β[1] = log(f0)
  nb_off = n * b
  grad = zeros(p)               # the score ∑ w (n_j − μ_j)(c_j − x)^l, l = 0 … degree
  momμ = zeros(2 * degree + 1)  # the μ-weighted design moments ∑ w μ_j (c_j − x)^k; the Gram is their Hankel
  for _ in 1:25
    fill!(grad, 0.0)
    fill!(momμ, 0.0)
    @inbounds for j in eachindex(centers)
      d = centers[j] - x
      w = _kernel_weight(centers[j], x, hx, kernel)
      w * hx < 1e-10 && continue        # cut on the raw kernel value (scale-free), not the 1/hx-folded weight
      η = β[1]
      dl = 1.0
      for l in 2:p
        dl *= d
        η += β[l] * dl
      end
      μ = nb_off * exp(η)
      wr = w * (counts[j] - μ)
      wμ = w * μ
      grad[1] += wr
      momμ[1] += wμ
      dk = 1.0
      for k in 1:(2 * degree)
        dk *= d
        momμ[k + 1] += wμ * dk
        k ≤ degree && (grad[k + 1] += wr * dk)
      end
    end
    δ = _hankel_solve(momμ, grad)
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
    density(summary::Histogram1D, x; h = nothing, kernel = gaussian_kernel, degree = 0,
            adaptive = false, support = nothing, start = nothing, normalize = false,
            binning = :linear)   # ≡ f_X

**Estimand.** The marginal (Lebesgue) density ``f_X = \\mathrm{d}\\mu_X/\\mathrm{d}\\lambda`` of the
conditioning variable ``X``, defined by ``P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda``.

**Estimator.** The smoothed estimate ``\\hat f_h(x)`` from the accumulated `summary` — a binned local
smoother of the counts ``n_j`` at bin locations ``\\ell_j`` (the within-bin centroid by default; see the
`binning` option), a KDE / WARPing, or a local-log-polynomial fit of ``\\log f`` at `degree` ``\\ge 1``.

# Arguments

- `summary` — the [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) of ``X``.
- `x` — the **query point(s)**: a scalar for the density at one point, or a vector for the whole curve.
- `h` — the **smoothing bandwidth**: the kernel scale, the one option that sets the bias–variance
  trade-off (larger `h` smooths more — less variance, more bias). Pass a positive number, or leave it
  `nothing` (default) to select one from the data by the kernel- and degree-aware **normal-reference
  rule** ([`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) — Silverman's
  ``(4/3)^{1/5}\\hat\\sigma n^{-1/5}`` for the Gaussian kernel at degree 0/1, wider at the ``n^{-1/9}``
  rate for degree 2). For a shape-adaptive choice, compute `bandwidth(summary; rule = :sj)` (or `:isj` /
  `:aic`) yourself and pass the number. The asymmetric [`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) / [`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) take their own scale
  ``b`` here (rate ``n^{-2/5}``) in place of the symmetric `h`.
- `kernel` — the smoothing weight ``K`` (any callable with ``\\int K = 1``), default the **symmetric**
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel). The compact-support symmetric
  [`epanechnikov_kernel`](@ref GBGStreamingStatistics.Kernels.epanechnikov_kernel) and, for **bounded /
  positive support**, the **asymmetric** `gamma_kernel` (``[0,\\infty)``, e.g. a droplet radius) and
  `beta_kernel` (``[0,1]``) — boundary-bias-free by
  construction — are catalogued with their supports and constants on [the Kernels page](@ref "Kernels");
  supply a custom [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel) for your own
  bounded support. An asymmetric kernel is its own degree-0, boundary-free estimator, so `degree`,
  `adaptive`, `support`, and a parametric `start` do not combine with it, and it takes its own scale ``b``
  in place of `h`.
- `degree` — the local-polynomial order. `degree = 0` (default) is the binned KDE / WARPing (weighted averaging of rounded points)
  ``\\hat f(x) = \\frac{1}{n} \\sum_j (n_j/h_j)  K((x-\\ell_j)/h_j)`` (``\\int K = 1``). `degree` ``\\ge 1`` is the
  **local-log-polynomial** (local-likelihood) estimator (Loader 1996): fit ``\\log f`` by a
  degree-`degree` polynomial through a binned **local Poisson** likelihood and report
  ``\\hat f(x) = e^{\\hat\\beta_0}``. `degree = 1` (log-linear) is **exact for an exponential density**,
  so it tracks the exponential / stretched-exponential **tails** a flat KDE smears; `degree = 2`
  (log-quadratic) is exact for a Gaussian. It falls back to the KDE where the local fit is degenerate or
  runs away in a sparse tail.
- `adaptive` — whether the bandwidth **varies with the local density**. `false` (default) holds `h`
  fixed everywhere; `true` uses the **Abramson** variable bandwidth — the kernel width scales by
  ``(\\tilde f/g)^{-1/2}`` (a pilot ``\\tilde f`` at the bandwidth `h`, ``g`` its geometric mean), widening
  in the sparse tails for far less tail variance (Abramson 1982). It is **orthogonal to how `h` is
  chosen**, so it composes with the data-driven bandwidth — pass `h = bandwidth(summary; rule = :sj)` to
  vary a *shape-adaptive* global scale locally (the chosen `h` becomes the Abramson pilot). (The degree-0
  KDE uses the sample-point ``h_j`` per bin; the local-log fit uses one balloon ``h(x)`` at the query.)
- `support` — the support of the density, passed as the tuple `(lower, upper)` (e.g. `(0.0, Inf)` for
  a positive quantity such as a droplet radius), or `nothing` (default, unbounded). When given, the
  degree-0 KDE is **renormalized** by the kernel mass inside the support: ``\\hat f(x) \\big/ \\int_a^b K_h(x-u)\\,\\mathrm{d}u`` (Jones 1993). This corrects the boundary bias that otherwise drops the estimate
  to ``\\approx`` half the true value at the edge. (`degree` ``\\ge 1`` is already boundary-adaptive, so
  `support` is a degree-0 option.)
- `start` — a **parametric start** (Hjort & Glad 1995), or `nothing` (default). Given one, `density`
  returns the semiparametric estimate
  ``\\hat f(x) = p(x;\\hat\\theta) \\frac{1}{nh}\\sum_j \\frac{n_j}{p(\\ell_j;\\hat\\theta)} K((x-\\ell_j)/h)`` —
  a moment-fitted model ``p(\\cdot;\\hat\\theta)`` times a KDE of its *correction*, **near-exact wherever
  the start is correct** (the correction is then ``\\equiv 1``) and degrading only to a bias-corrected KDE
  where it is not. ``\\hat\\theta`` is fit by method of moments from the recovered ``\\hat\\mu``,
  ``\\hat\\sigma^2``, so the start is itself a functional of the summary. Pass a built-in `normal_start`,
  `lognormal_start`, `gamma_start`, `exponential_start`, or `inverse_gaussian_start` — [the Parametric
  starts section](@ref "Parametric starts") lists each with its moment fit — or supply your own
  [`ParametricStart`](@ref GBGStreamingStatistics.MarginalDensity.ParametricStart) (a moment fit, a
  log-density, a support). It is a degree-0, fixed-`h`, symmetric-kernel estimator (`degree` ``\\ge 1``,
  `adaptive`, `support`, and the asymmetric kernels do not combine with it; use `normalize` for
  ``\\int\\hat f = 1``).
- `normalize` — when `true`, divide the estimate by ``\\int\\hat f`` over the summary's range
  (intersected with the estimator's support), so ``\\int\\hat f = 1`` **exactly**. The symmetric KDE
  already integrates to ``\\approx 1``; the asymmetric kernels and a parametric start do not, so this is
  their finishing step.
- `binning` — where each bin's mass sits for the post-hoc evaluation. `:linear` (the default) places
  it at the recovered within-bin **centroid** ``\\ell_j = c_j + S_1^{(j)}/n_j`` (from the stored first power-sum
  ``S_1``), which matches the bin's own first moment, so the binned estimate converges to the exact KDE at
  ``O(\\Delta^2)`` in the bin width ``\\Delta``; `:simple` places it at the geometric bin center ``\\ell_j = c_j``,
  the ``O(\\Delta)`` mass-at-center approximation. The centroid is read from data already in the summary, so
  the higher-order accuracy is free; it applies to every path (KDE, local-log, start, asymmetric kernel).

The four corners (fixed/variable ``\\times`` KDE/local-log) and everything between are one call. A vector
of query points returns a vector.

```jldoctest
julia> h = Histogram1D(-1.0, 1.0, 2);

julia> observe!(h, [0.4]);          # one sample at 0.4

julia> round(density(h, 0.4; h = 0.2); digits = 4)   # linear binning centers the bump at 0.4 (the centroid)
1.9947
```

# References

- I. S. Abramson (1982), [*On bandwidth variation in kernel estimates—a square root law*](https://doi.org/10.1214/aos/1176345986), Ann. Statist. 10(4) (variable bandwidth).
- C. R. Loader (1996), [*Local likelihood density estimation*](https://doi.org/10.1214/aos/1032298287), Ann. Statist. 24(4) (the local-log-polynomial estimator).
- M. C. Jones (1993), [*Simple boundary correction for kernel density estimation*](https://doi.org/10.1007/BF00147776), Statist. Comput. 3(3) (the `support` renormalization).
- N. L. Hjort & I. K. Glad (1995), [*Nonparametric density estimation with a parametric start*](https://doi.org/10.1214/aos/1176324627), Ann. Statist. 23(3) (the `start` estimator).
"""
function density(summary::Histogram1D, x::Real; h = nothing, kernel = gaussian_kernel,
                 degree::Integer = 0, adaptive::Bool = false, support = nothing,
                 start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  if normalize
    Z = _normalization_constant(summary; degree = degree, h = h, adaptive = adaptive,
                                kernel = kernel, support = support, start = start, binning = binning)
    fx = density(summary, x; degree = degree, h = h, adaptive = adaptive, kernel = kernel,
                 support = support, start = start, normalize = false, binning = binning)
    return Z > 0 ? fx / Z : fx
  end
  if start !== nothing
    _check_start(kernel, degree, adaptive, support)
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
    (degree == 0 && !adaptive && support === nothing) ||
      throw(ArgumentError("the asymmetric kernel is its own degree-0, boundary-free estimator; " *
                          "`degree`, `adaptive`, and `support` do not apply"))
    n = _total(summary)
    iszero(n) && return 0.0
    bb = h === nothing ? _asymmetric_default_b(summary) : float(h)
    bb > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return _kde_h(counts, centers, n, x, bb, kernel)   # one degree-0 evaluator, any kernel protocol
  end
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  (support === nothing || degree == 0) ||
    throw(ArgumentError("`support` (boundary renormalization) applies to the degree-0 KDE; " *
                        "degree ≥ 1 (local-log) is already boundary-adaptive (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return 0.0
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = _bin_locations(summary, counts, binning)
  if degree == 0 && !adaptive
    fx = _kde_h(counts, centers, n, x, hh, kernel)
    return support === nothing ? fx : fx / _kernel_mass(kernel, x, hh, support[1], support[2])
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    fx = _kde_adaptive(counts, centers, n, x, hj, kernel)
    return support === nothing ? fx : fx / _kernel_mass(kernel, x, hh, support[1], support[2])
  end
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  f0 = _kde_h(counts, centers, n, x, hh, kernel)
  return _local_log(counts, centers, n, x, _balloon(hh, g, f0, adaptive), kernel, degree,
                    inv(summary.scale), f0)
end

function density(summary::Histogram1D, xs::AbstractVector; h = nothing,
                 kernel = gaussian_kernel, degree::Integer = 0, adaptive::Bool = false,
                 support = nothing, start = nothing, normalize::Bool = false,
                 binning::Symbol = :linear)
  if normalize
    Z = _normalization_constant(summary; degree = degree, h = h, adaptive = adaptive,
                                kernel = kernel, support = support, start = start, binning = binning)
    fxs = density(summary, xs; degree = degree, h = h, adaptive = adaptive, kernel = kernel,
                  support = support, start = start, normalize = false, binning = binning)
    return Z > 0 ? fxs ./ Z : fxs
  end
  if start !== nothing
    _check_start(kernel, degree, adaptive, support)
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
    (degree == 0 && !adaptive && support === nothing) ||
      throw(ArgumentError("the asymmetric kernel is its own degree-0, boundary-free estimator; " *
                          "`degree`, `adaptive`, and `support` do not apply"))
    n = _total(summary)
    iszero(n) && return zeros(Float64, length(xs))
    bb = h === nothing ? _asymmetric_default_b(summary) : float(h)
    bb > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    return [_kde_h(counts, centers, n, x, bb, kernel) for x in xs]
  end
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  (support === nothing || degree == 0) ||
    throw(ArgumentError("`support` (boundary renormalization) applies to the degree-0 KDE; " *
                        "degree ≥ 1 (local-log) is already boundary-adaptive (got degree = $degree)"))
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = _bin_locations(summary, counts, binning)
  if degree == 0 && !adaptive
    fxs = [_kde_h(counts, centers, n, x, hh, kernel) for x in xs]
    return support === nothing ? fxs : [f / _kernel_mass(kernel, x, hh, support[1], support[2]) for (f, x) in zip(fxs, xs)]
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    fxs = [_kde_adaptive(counts, centers, n, x, hj, kernel) for x in xs]
    return support === nothing ? fxs : [f / _kernel_mass(kernel, x, hh, support[1], support[2]) for (f, x) in zip(fxs, xs)]
  end
  b = inv(summary.scale)
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  return [begin
            f0 = _kde_h(counts, centers, n, x, hh, kernel)
            _local_log(counts, centers, n, x, _balloon(hh, g, f0, adaptive), kernel, degree, b, f0)
          end for x in xs]
end


# --- pointwise standard error (§4.1): SE(x) = √(f̂(x)) · scale, the KDE variability band ---
# The boundary-restricted roughness R_a(x) = ∫_{[a,b]} K((x−u)/h)²/h du (the ∫K² analog of `_kernel_mass`'s
# a(x)): the CDF of K² up to (x−a)/h minus up to (x−b)/h. Interior / no `support` ⇒ the full R(K) = ∫K².
_kernel_sq_cdf(kernel, t) = _quad(s -> s ≤ t ? kernel(s)^2 : 0.0)

function _kernel_sq_mass(kernel, x, h, a, b)
  full = _quad(s -> kernel(s)^2)
  hi = isfinite(a) ? _kernel_sq_cdf(kernel, (x - a) / h) : full
  lo = isfinite(b) ? _kernel_sq_cdf(kernel, (x - b) / h) : 0.0
  return hi - lo
end


# The variability band is for a symmetric kernel at degree 0/1 (the asymmetric-kernel variance and
# degree ≥ 2 are recorded coverage gaps).
function _check_se(kernel, degree)
  kernel isa AsymmetricKernel &&
    throw(ArgumentError("the density variability band is for a symmetric kernel; it is not supported for the asymmetric gamma/beta kernels"))
  (degree == 0 || degree == 1) ||
    throw(ArgumentError("the density variability band is degree-0/1 (KDE); it is not supported for degree ≥ 2 (got degree = $degree)"))
end


# The f̂-independent SE scale, so SE(x) = √(f̂(x)) · scale = √(R_eff(x)·f̂(x)/(n·h_eff(x)))/a(x):
#   h_eff  the Abramson balloon bandwidth h(x) when `adaptive` (else h), so the variance uses the *local*
#          bandwidth in the tail; R_eff, a(x)  the boundary-restricted roughness / kernel mass when `support`
#          (else R(K), 1), the Jones renormalization's variance inflation near an edge. `start` leaves the
#          leading variance unchanged (θ̂'s O(1/n) uncertainty is a documented omission), so needs nothing.
function _se_scale(summary::Histogram1D, x, hh, kernel, adaptive, support, binning)
  n = _total(summary)
  heff = hh
  if adaptive
    counts = Array(summary.counts)
    centers = _bin_locations(summary, counts, binning)
    _, g = _abramson(counts, centers, n, hh, kernel)
    pilot = _kde_h(counts, centers, n, x, hh, kernel)
    heff = _balloon(hh, g, pilot, true)
  end
  slower, supper = support === nothing ? (-Inf, Inf) : (float(support[1]), float(support[2]))
  Reff = _kernel_sq_mass(kernel, x, hh, slower, supper)
  a = _kernel_mass(kernel, x, hh, slower, supper)
  return sqrt(Reff / (n * heff)) / a
end


# SE(x) = √(max(Var, 0)), the finite-sample KDE variance Var = R_eff·f̂/(n·h_eff·a²) − f̂²/n: the leading
# (nh)⁻¹ term (`_se_scale`²) minus the n⁻¹ term (Wand & Jones 1995) that keeps the band from over-widening
# where f̂ is not small (the mode, the boundary), given the recovered f̂ = `fx`.
function _se_from(summary::Histogram1D, x, fx, hh, kernel, adaptive, support, binning)
  n = _total(summary)
  s2 = _se_scale(summary, x, hh, kernel, adaptive, support, binning)^2   # R_eff / (n·h_eff·a²)
  f = max(fx, 0.0)
  return sqrt(max(s2 * f - f^2 / n, 0.0))
end


"""
    density_standard_error(summary::Histogram1D, x;
                           h = nothing, kernel = gaussian_kernel, degree = 0,
                           adaptive = false, support = nothing, start = nothing,
                           normalize = false, binning = :linear)   # ≡ SE[f̂(x)]

**Estimand.** The pointwise standard error of the density estimator
`density` at ``x`` — ``\\mathrm{SE}(x) =
\\sqrt{\\operatorname{Var}\\hat f_h(x)}`` with the kernel-density asymptotic variance

```math
\\operatorname{Var}\\hat f_h(x) \\approx \\frac{R(K) f_X(x)}{n h} - \\frac{f_X(x)^2}{n}, \\qquad R(K) = \\int K^2.
```

**Estimator.** ``\\sqrt{R(K) \\hat f_h(x)/(n h) - \\hat f_h(x)^2/n}`` from the recovered ``\\hat f_h(x)``,
the count ``n``, the bandwidth `h` (default the normal reference), and the roughness ``R(K) = \\int K^2`` of the
smoothing `kernel` ``K`` — any symmetric ``\\int K = 1`` callable, default
[`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel) — by quadrature. The
pointwise **variability band** is ``\\hat f_h(x) \\pm z \\mathrm{SE}(x)`` (``z = 1.96`` for 95\\%) — a
band about ``\\mathbb{E}\\hat f_h``, so it covers ``f_X`` where the ``O(h^2)`` smoothing bias is small
and narrows that coverage at a sharp mode or in a tail.

**It tracks the estimator's options** — pass the *same* options you gave `density` and the band is the
variance of *that* estimate (or use
[`density_with_standard_error`](@ref
GBGStreamingStatistics.MarginalDensity.density_with_standard_error) to get the estimate and its band from
one call). The variance is adjusted per option: `adaptive` replaces ``h`` by the **local** Abramson
bandwidth ``h(x)`` (a wider window — hence smaller variance — in the tail; the sample-point pilot's own
variability is an unmodelled higher-order term, so the band runs slightly narrow in the deep tail); `support` (Jones
renormalization) replaces ``R(K)`` by the boundary-restricted roughness ``R_a(x) =
\\frac{1}{h}\\int_{\\text{support}} K((x-u)/h)^2\\,\\mathrm du`` and divides by the squared in-support kernel
mass ``a(x)^2``, inflating the variance near an edge; a parametric `start` leaves the leading variance
``R(K)\\hat f/(nh)`` unchanged (the fitted-parameter uncertainty is an ``O(1/n)`` omission). Restricted to
`degree` 0/1 and a **symmetric** kernel (``\\texttt{degree} \\ge 2`` and an asymmetric gamma/beta kernel
are not supported); a vector of query points returns a vector.

Arguments: as for [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) — the options specify the
estimator whose standard error is returned.

```jldoctest
julia> h = Histogram1D(-6.0, 6.0, 600);

julia> observe!(h, [-0.4, 0.1, 0.2, 0.9, 1.3]);

julia> density_standard_error(h, 0.5) > 0
true
```

Reference: M. P. Wand & M. C. Jones (1995), [*Kernel Smoothing*](https://doi.org/10.1201/b14876), Chapman & Hall.
"""
function density_standard_error(summary::Histogram1D, x::Real; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, adaptive::Bool = false,
                    support = nothing, start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  _check_se(kernel, degree)
  n = _total(summary)
  iszero(n) && return 0.0
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  fx = density(summary, x; h = hh, kernel = kernel, degree = degree, adaptive = adaptive,
               support = support, start = start, normalize = normalize, binning = binning)
  return _se_from(summary, x, fx, hh, kernel, adaptive, support, binning)
end

function density_standard_error(summary::Histogram1D, xs::AbstractVector; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, adaptive::Bool = false,
                    support = nothing, start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  _check_se(kernel, degree)
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  fxs = density(summary, xs; h = hh, kernel = kernel, degree = degree, adaptive = adaptive,
                support = support, start = start, normalize = normalize, binning = binning)
  return [_se_from(summary, x, f, hh, kernel, adaptive, support, binning)
          for (f, x) in zip(fxs, xs)]
end


"""
    density_with_standard_error(summary::Histogram1D, x;
                                h = nothing, kernel = gaussian_kernel, degree = 0,
                                adaptive = false, support = nothing, start = nothing,
                                normalize = false, binning = :linear)

The `density` estimate **and** its pointwise
[`density_standard_error`](@ref GBGStreamingStatistics.MarginalDensity.density_standard_error), from a
**single option specification** — so the estimate and its error describe the same estimator by
construction (the estimate is computed once and reused). Returns the `NamedTuple`
`(value = f̂(x), standard_error = SE(x))`; a vector `x` gives a `(value, standard_error)` of vectors. The
standard error's estimand and per-option variance are documented at `density_standard_error`.

Arguments: as for [`density`](@ref GBGStreamingStatistics.MarginalDensity.density).

```jldoctest
julia> h = Histogram1D(-6.0, 6.0, 600);

julia> observe!(h, [-0.4, 0.1, 0.2, 0.9, 1.3]);

julia> r = density_with_standard_error(h, 0.5); (r.value > 0, r.standard_error > 0)
(true, true)
```
"""
function density_with_standard_error(summary::Histogram1D, x::Real; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, adaptive::Bool = false,
                    support = nothing, start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  _check_se(kernel, degree)
  n = _total(summary)
  iszero(n) && return (value = 0.0, standard_error = 0.0)
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  fx = density(summary, x; h = hh, kernel = kernel, degree = degree, adaptive = adaptive,
               support = support, start = start, normalize = normalize, binning = binning)
  return (value = fx,
          standard_error = _se_from(summary, x, fx, hh, kernel, adaptive, support, binning))
end

function density_with_standard_error(summary::Histogram1D, xs::AbstractVector; h = nothing,
                    kernel = gaussian_kernel, degree::Integer = 0, adaptive::Bool = false,
                    support = nothing, start = nothing, normalize::Bool = false, binning::Symbol = :linear)
  _check_se(kernel, degree)
  n = _total(summary)
  iszero(n) && return (value = zeros(Float64, length(xs)), standard_error = zeros(Float64, length(xs)))
  hh = h === nothing ? _normal_reference_h(summary, kernel, degree) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  fxs = density(summary, xs; h = hh, kernel = kernel, degree = degree, adaptive = adaptive,
                support = support, start = start, normalize = normalize, binning = binning)
  ses = [_se_from(summary, x, f, hh, kernel, adaptive, support, binning)
         for (f, x) in zip(fxs, xs)]
  return (value = fxs, standard_error = ses)
end

end # module MarginalDensity
