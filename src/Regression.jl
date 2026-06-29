# Regression — §4.2: the conditional mean m(x) and conditional variance σ²(x), by local-linear
# regression from the joint Histogram2D. The estimand-first reference module — write every estimator
# docstring to this pattern (`docs/ARCHITECTURE.md` §2, §4.2).

module Regression

using LinearAlgebra

using ..Histograms: Histogram2D, _xbincenter, _ybincenter

export conditional_mean, conditional_variance


# Unnormalised Gaussian smoothing weight; the 1/√(2π) and the 1/h both cancel in the local-linear
# ratio below, so they are dropped.
@inline _gaussian(u) = exp(-u^2 / 2)


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


# Local-linear fit at `x`: the kernel-weighted sums ``T_\ell = \sum_{i_x} K_h(c^x-x)(c^x-x)^\ell
# N_{i_x}`` and ``P_\ell = \sum_{i_x} K_h(c^x-x)(c^x-x)^\ell R_{i_x}`` (R the per-column response
# total), solved as ``(T_2 P_0 - T_1 P_1)/(T_0 T_2 - T_1^2)``. Returns `NaN` where there is no local
# data (degenerate design).
function _local_linear(cx, N, R, x, invh)
  T0 = 0.0; T1 = 0.0; T2 = 0.0; P0 = 0.0; P1 = 0.0
  @inbounds for ix in eachindex(cx)
    d = cx[ix] - x
    w = _gaussian(d * invh)
    wN = w * N[ix]
    wR = w * R[ix]
    T0 += wN; T1 += wN * d; T2 += wN * d * d
    P0 += wR; P1 += wR * d
  end
  denom = T0 * T2 - T1 * T1
  return denom > 0 ? (T2 * P0 - T1 * P1) / denom : oftype(denom, NaN)
end


"""
    conditional_mean(summary::Histogram2D, x; h)   # ≡ m(x)

**Estimand.** the conditional expectation
``m(x) = \\mathbb{E}[Y \\mid X = x] = \\int y\\, f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y`` — the
``L^2(\\sigma(X))`` projection of the response ``Y`` onto the conditioning variable ``X``, i.e. the
regression function.

**Estimator** (local-linear; Fan & Gijbels). From the joint `summary` (a
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X, Y)``) form the
kernel-weighted sums ``T_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell`` and
``P_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell Y_i`` — recovered from the per-X-column counts and
response totals ``\\sum Y`` (themselves exact, from the per-cell power-sums) — and solve the local
normal equations:

```math
\\hat m_h(x) = \\frac{T_2 P_0 - T_1 P_1}{T_0 T_2 - T_1^{2}}.
```

Local-linear is the default — it removes the design-density bias ``2 m'(x) f'(x)/f(x)`` of the
Nadaraya–Watson reduction ``\\hat m_h = P_0/T_0``. The bandwidth ``h`` is supplied (the response-aware
local plug-in ``h^\\star(x)`` is forthcoming, §4.5); a vector of query points returns a vector.
Consistency (`docs/FORMALIZATION.md`): ``\\operatorname{bias}\\hat m_h(x) \\approx \\tfrac12 h^2
\\mu_2(K)\\, m''(x)``, ``\\operatorname{Var}\\hat m_h(x) \\approx R(K)\\,\\sigma^2(x)/(n\\,h\\,f_X(x))``.

Reference: J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its Applications*](https://doi.org/10.1201/9780203748725),
Chapman & Hall.
"""
function conditional_mean(summary::Histogram2D, x::Real; h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  cx, N, sumY, _ = _column_sums(summary)
  return _local_linear(cx, N, sumY, x, inv(float(h)))
end

function conditional_mean(summary::Histogram2D, xs::AbstractVector; h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  cx, N, sumY, _ = _column_sums(summary)
  invh = inv(float(h))
  return [_local_linear(cx, N, sumY, x, invh) for x in xs]
end


"""
    conditional_variance(summary::Histogram2D, x; h)   # ≡ σ²(x)

**Estimand.** the conditional variance
``\\sigma^2(x) = \\operatorname{Var}(Y \\mid X = x) = \\mathbb{E}[Y^2 \\mid X = x] - m(x)^2``, with
``m(x)`` the [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean).

**Estimator.** the same local-linear fit applied to the response ``Y^2`` (its per-X-column totals
``\\sum Y^2``, exact from the per-cell power-sums), as the plug-in ``\\widehat{\\mathbb{E}}[Y^2\\mid x] -
\\hat m_h(x)^2`` clamped to be non-negative. This simple difference carries the local-linear smoothing
bias of the curved ``\\mathbb{E}[Y^2\\mid x]``; the bias-reducing **residual-based** estimator (a second
local fit to the squared residuals; Fan & Yao 1998) and a one-step robust local M-estimator for
heavy-tailed channels (e.g. ``(\\Delta X)^6``; Fan & Jiang 2000) are forthcoming.

References:

- J. Fan & Q. Yao (1998), [*Efficient estimation of conditional variance functions in stochastic regression*](https://doi.org/10.1093/biomet/85.3.645), Biometrika 85(3) (the residual-based estimator).
- J. Fan & J. Jiang (2000), [*Variable bandwidth and one-step local M-estimator*](https://doi.org/10.1007/BF02903849), Sci. China Ser. A 43 (the robust refinement).
"""
function conditional_variance(summary::Histogram2D, x::Real; h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  cx, N, sumY, sumY2 = _column_sums(summary)
  invh = inv(float(h))
  m = _local_linear(cx, N, sumY, x, invh)
  ey2 = _local_linear(cx, N, sumY2, x, invh)
  v = ey2 - m^2
  return v > 0 ? v : oftype(v, isnan(v) ? NaN : 0.0)
end

function conditional_variance(summary::Histogram2D, xs::AbstractVector; h::Real)
  return [conditional_variance(summary, x; h = h) for x in xs]
end

end # module Regression
