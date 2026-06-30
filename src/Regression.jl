# Regression — §4.2: the conditional mean m(x) and conditional variance σ²(x), by local-polynomial
# regression from the joint Histogram2D. The estimand-first reference module — write every estimator
# docstring to this pattern (`docs/ARCHITECTURE.md` §2, §4.2). The kernel-weighted local-polynomial
# fit is the shared `LocalPolynomial.local_fit` primitive (§11); the conditional mean is its β₀.

module Regression

using ..Histograms: Histogram2D, _xbincenter, _ybincenter
using ..LocalPolynomial: local_fit
using ..DensityEstimation: gaussian_kernel

export conditional_mean, conditional_variance


# Local-polynomial order admitted by the conditional-moment estimators: 0 (Nadaraya–Watson),
# 1 (local-linear, the default), 2 (local-quadratic — exposes the curvature m'' that h⋆(x) needs).
@inline function _check_degree(degree)
  (degree == 0 || degree == 1 || degree == 2) ||
    throw(ArgumentError("degree must be 0, 1, or 2 (got degree = $degree)"))
end


# Per-X-column reductions over a Histogram2D: for each X-bin, the count ``N``, the response total
# ``\sum Y`` and ``\sum Y^2`` — recovered *exactly* from the per-cell counts and power-sums (``S_{01}
# = \sum (Y-c^y)`` and ``S_{02}`` use the actual ``Y``, not the cell center). Host-side; a no-op
# `Array` for CPU summaries.
function _column_sums(summary::Histogram2D)
  counts = Array(summary.counts)
  S = Array(summary.powersums)
  nx = summary.nx
  ny = summary.ny
  cx = [_xbincenter(summary, ix) for ix in 1:nx]
  N = zeros(Float64, nx)
  sumY = zeros(Float64, nx)
  sumY2 = zeros(Float64, nx)
  for iy in 1:ny
    cy = _ybincenter(summary, iy)
    for ix in 1:nx
      cell = ix + (iy - 1) * nx
      nc = counts[cell]
      iszero(nc) && continue
      s01 = S[2, cell]
      s02 = S[5, cell]
      N[ix] += nc
      sumY[ix] += cy * nc + s01                       # Σ Y over the column, exact
      sumY2[ix] += cy * cy * nc + 2 * cy * s01 + s02  # Σ Y² over the column, exact
    end
  end
  return cx, N, sumY, sumY2
end


"""
    conditional_mean(summary::Histogram2D, x;
                     h, kernel = gaussian_kernel, degree = 1)   # ≡ m(x)

**Estimand.** the conditional expectation
``m(x) = \\mathbb{E}[Y \\mid X = x] = \\int y\\, f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y`` — the
``L^2(\\sigma(X))`` projection of the response ``Y`` onto the conditioning variable ``X``, i.e. the
regression function.

**Estimator** (local-polynomial; Fan & Gijbels). From the joint `summary` (a
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X, Y)``) form the
kernel-weighted design and response sums ``T_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell`` and
``P_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell Y_i`` — recovered from the per-X-column counts and
response totals ``\\sum Y`` (themselves exact, from the per-cell power-sums) — and solve the local
normal equations of order `degree` for the jet; the estimate is the fitted value
``\\hat m_h(x) = \\beta_0``. At the default `degree = 1` (local-linear) this is

```math
\\hat m_h(x) = \\frac{T_2 P_0 - T_1 P_1}{T_0 T_2 - T_1^{2}},
```

which removes the design-density bias ``2 m'(x) f'(x)/f(x)`` of the Nadaraya–Watson reduction
``\\hat m_h = P_0/T_0`` (the `degree = 0` local-constant fit). The bandwidth ``h`` is supplied (the
response-aware local plug-in ``h^\\star(x)`` is forthcoming, §4.5); `degree = 2` (local-quadratic)
additionally exposes the curvature ``m''(x)`` that ``h^\\star(x)`` needs. A vector of query points
returns a vector. Consistency (`docs/FORMALIZATION.md`): ``\\operatorname{bias}\\hat m_h(x) \\approx
\\tfrac12 h^2 \\mu_2(K)\\, m''(x)``, ``\\operatorname{Var}\\hat m_h(x) \\approx
R(K)\\,\\sigma^2(x)/(n\\,h\\,f_X(x))``.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``m`` is evaluated (a scalar for one point, or a vector
  for the whole curve).
- `h` — the smoothing bandwidth (``h > 0``), supplied.
- `kernel` — the smoothing kernel ``K`` (any callable with ``\\int K = 1``; default `gaussian_kernel`).
- `degree` — the local-polynomial order: `0` (Nadaraya–Watson), `1` (local-linear, the default), or
  `2` (local-quadratic).

Reference: J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its Applications*](https://doi.org/10.1201/9780203748725),
Chapman & Hall.
"""
function conditional_mean(summary::Histogram2D, x::Real; h::Real,
                          kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumY, _ = _column_sums(summary)
  return local_fit(cx, N, sumY, x, h, kernel, degree)[1]
end

function conditional_mean(summary::Histogram2D, xs::AbstractVector; h::Real,
                          kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumY, _ = _column_sums(summary)
  return [local_fit(cx, N, sumY, x, h, kernel, degree)[1] for x in xs]
end


"""
    conditional_variance(summary::Histogram2D, x;
                         h, kernel = gaussian_kernel, degree = 1)   # ≡ σ²(x)

**Estimand.** the conditional variance
``\\sigma^2(x) = \\operatorname{Var}(Y \\mid X = x) = \\mathbb{E}[Y^2 \\mid X = x] - m(x)^2``, with
``m(x)`` the [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean).

**Estimator.** the same local-polynomial fit applied to the response ``Y^2`` (its per-X-column totals
``\\sum Y^2``, exact from the per-cell power-sums), as the plug-in ``\\widehat{\\mathbb{E}}[Y^2\\mid x] -
\\hat m_h(x)^2`` clamped to be non-negative. This simple difference carries the local-polynomial
smoothing bias of the curved ``\\mathbb{E}[Y^2\\mid x]``; the bias-reducing **residual-based** estimator
(a second local fit to the squared residuals; Fan & Yao 1998) and a one-step robust local M-estimator
for heavy-tailed channels (e.g. ``(\\Delta X)^6``; Fan & Jiang 2000) are forthcoming.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``\\sigma^2`` is evaluated (scalar or vector).
- `h` — the smoothing bandwidth (``h > 0``), supplied.
- `kernel` — the smoothing kernel ``K`` (any callable with ``\\int K = 1``; default `gaussian_kernel`).
- `degree` — the local-polynomial order for the two underlying fits: `0`, `1` (the default), or `2`.

References:

- J. Fan & Q. Yao (1998), [*Efficient estimation of conditional variance functions in stochastic regression*](https://doi.org/10.1093/biomet/85.3.645), Biometrika 85(3) (the residual-based estimator).
- J. Fan & J. Jiang (2000), [*Variable bandwidth and one-step local M-estimator*](https://doi.org/10.1007/BF02903849), Sci. China Ser. A 43 (the robust refinement).
"""
function conditional_variance(summary::Histogram2D, x::Real; h::Real,
                              kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumY, sumY2 = _column_sums(summary)
  m = local_fit(cx, N, sumY, x, h, kernel, degree)[1]
  ey2 = local_fit(cx, N, sumY2, x, h, kernel, degree)[1]
  v = ey2 - m^2
  return v > 0 ? v : oftype(v, isnan(v) ? NaN : 0.0)
end

function conditional_variance(summary::Histogram2D, xs::AbstractVector; h::Real,
                              kernel = gaussian_kernel, degree::Integer = 1)
  return [conditional_variance(summary, x; h = h, kernel = kernel, degree = degree) for x in xs]
end

end # module Regression
