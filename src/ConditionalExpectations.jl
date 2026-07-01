# ConditionalExpectations — §4.2: the conditional mean m(x) and conditional variance σ²(x), by local-polynomial
# regression from the joint Histogram2D. The estimand-first reference module — write every estimator
# docstring to this pattern (`docs/ARCHITECTURE.md` §2, §4.2). The kernel-weighted local-polynomial
# fit is the shared `LocalPolynomial.local_fit` primitive (§11); the conditional mean is its β₀.

module ConditionalExpectations

using ..Histograms: Histogram2D, _column_sums
using ..LocalPolynomial: local_fit
using ..Kernels: gaussian_kernel

export conditional_mean, conditional_variance


# Local-polynomial order admitted by the conditional-moment estimators: 0 (Nadaraya–Watson),
# 1 (local-linear, the default), 2 (local-quadratic — exposes the curvature m'' that h⋆(x) needs).
@inline function _check_degree(degree)
  (degree == 0 || degree == 1 || degree == 2) ||
    throw(ArgumentError("degree must be 0, 1, or 2 (got degree = $degree)"))
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
``\\hat m_h = P_0/T_0`` (the `degree = 0` local-constant fit). The bandwidth ``h`` is supplied (or the
response-aware local plug-in ``h^\\star(x)`` from `bandwidth(summary, x; rule = :local)`);
`degree = 2` (local-quadratic) additionally exposes the curvature ``m''(x)`` that ``h^\\star(x)`` needs. A vector of query points
returns a vector. Consistency: ``\\operatorname{bias}\\hat m_h(x) \\approx
\\tfrac12 h^2 \\mu_2(K)\\, m''(x)``, ``\\operatorname{Var}\\hat m_h(x) \\approx
R(K)\\,\\sigma^2(x)/(n\\,h\\,f_X(x))``.

**Boundary adaptivity.** The local-linear fit (`degree = 1`, the default) is *boundary-adaptive* (Fan
1992): conditioning on a bounded or positive quantity — a droplet radius ``X \\ge 0``, say — is
first-order correct near the edge with the symmetric `kernel`, because the local-linear fit corrects the
``O(h)`` boundary bias to the interior ``O(h^2)`` automatically. The `degree = 0` Nadaraya–Watson fit does
**not** — it keeps the ``O(h)`` edge bias — so near a bounded edge **use `degree = 1`** (the default; being
local-linear it is also unbiased for a linear mean). A `degree = 0` asymmetric `kernel` (`gamma_kernel` /
`beta_kernel`) is *admitted* too — the uniform degree-0 boundary kernel, no jet, so `degree = 0` only — but
the Chen kernels are a **density** boundary tool whose asymmetry interacts with the regression design
density, so they do **not** generally improve the conditional mean; prefer `degree = 1` for a bounded
conditioning variable. A bounded *response* is the other case: there the response kernel can be asymmetric —
see [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density).

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``m`` is evaluated (a scalar for one point, or a vector
  for the whole curve).
- `h` — the smoothing bandwidth (``h > 0``), supplied.
- `kernel` — the smoothing kernel: a symmetric callable ``K`` (``\\int K = 1``; default `gaussian_kernel`,
  at any `degree`), or an asymmetric `gamma_kernel` / `beta_kernel` for a bounded/positive conditioning
  variable (`degree = 0` only).
- `degree` — the local-polynomial order: `0` (Nadaraya–Watson), `1` (local-linear, the default), or
  `2` (local-quadratic); an asymmetric `kernel` restricts it to `0`.

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


# Clamp a variance estimate to be non-negative; preserve NaN (off-support / degenerate design).
@inline _clamp_var(v) = v > 0 ? v : oftype(v, isnan(v) ? NaN : 0.0)


# Fan & Yao (1998) residual response for the bias-reduced conditional variance: per X-column, the total
# of squared residuals Σ(Y − m̂(X))², with m̂(X) approximated by the column-center fit μ_j = m̂(c^x_j)
# (exact to O(Δ), the bin width). The column total is then ΣY² − 2 μ_j ΣY + μ_j² N from the stored column
# sums. A second local fit of this response (≡ Ê[(Y−m̂)²|x]) estimates σ²(x) with the curved-E[Y²|x] bias
# of the plug-in removed — E[(Y−m)²|x] = σ²(x) is far flatter than E[Y²|x].
function _residual_response(cx, N, sumY, sumY2, h, kernel, degree)
  SR2 = similar(sumY2)
  @inbounds for j in eachindex(cx)
    if N[j] > 0
      μj = local_fit(cx, N, sumY, cx[j], h, kernel, degree)[1]
      SR2[j] = sumY2[j] - 2 * μj * sumY[j] + μj^2 * N[j]
    else
      SR2[j] = zero(eltype(SR2))
    end
  end
  return SR2
end


# σ²(x) at a single point by the chosen method, from the pre-reduced column sums.
function _conditional_variance(cx, N, sumY, sumY2, x, h, kernel, degree, method)
  if method === :residual
    SR2 = _residual_response(cx, N, sumY, sumY2, h, kernel, degree)
    return _clamp_var(local_fit(cx, N, SR2, x, h, kernel, degree)[1])
  elseif method === :plugin
    m = local_fit(cx, N, sumY, x, h, kernel, degree)[1]
    ey2 = local_fit(cx, N, sumY2, x, h, kernel, degree)[1]
    return _clamp_var(ey2 - m^2)
  end
  throw(ArgumentError("conditional_variance `method` must be :residual (Fan–Yao, default) or " *
                      ":plugin (got :$method)"))
end


"""
    conditional_variance(summary::Histogram2D, x;
                         h, kernel = gaussian_kernel,
                         degree = 1, method = :residual)   # ≡ σ²(x)

**Estimand.** the conditional variance
``\\sigma^2(x) = \\operatorname{Var}(Y \\mid X = x) = \\mathbb{E}[Y^2 \\mid X = x] - m(x)^2``, with
``m(x)`` the [`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean).

**Estimator** (`method = :residual`, the default — Fan & Yao 1998). The **two-stage residual-based**
estimator: fit ``\\hat m`` at each X-column center, then take a second local fit of the per-column
squared residuals, ``\\hat\\sigma^2(x) = \\widehat{\\mathbb{E}}\\bigl[(Y - \\hat m(X))^2 \\mid x\\bigr]``.
Because ``\\mathbb{E}[(Y - m(X))^2 \\mid x] = \\sigma^2(x)`` is far flatter than
``\\mathbb{E}[Y^2 \\mid x]``, its local fit is much less biased — it removes the local-polynomial
smoothing bias of the curved ``\\mathbb{E}[Y^2 \\mid x]`` that the plug-in carries. The per-column
squared-residual total is recovered exactly from the stored column sums,
``\\sum Y^2 - 2\\hat m \\sum Y + \\hat m^2 N``.

**Estimator** (`method = :plugin`). The simple difference ``\\widehat{\\mathbb{E}}[Y^2 \\mid x] -
\\hat m_h(x)^2`` (the same local fit applied to the ``Y^2`` response), clamped non-negative — cheaper
(no first stage), but it carries the curved-``\\mathbb{E}[Y^2 \\mid x]`` smoothing bias. A one-step
robust local M-estimator for heavy-tailed channels (e.g. ``(\\Delta X)^6``; Fan & Jiang 2000) is
forthcoming.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``\\sigma^2`` is evaluated (scalar or vector).
- `h` — the smoothing bandwidth (``h > 0``), supplied.
- `kernel` — the smoothing kernel ``K`` (any callable with ``\\int K = 1``; default `gaussian_kernel`).
- `degree` — the local-polynomial order for the underlying fits: `0`, `1` (the default), or `2`.
- `method` — `:residual` (the default, the bias-reduced Fan & Yao two-stage estimator) or `:plugin`
  (the simple ``\\widehat{\\mathbb{E}}[Y^2 \\mid x] - \\hat m^2`` difference).

References:

- J. Fan & Q. Yao (1998), [*Efficient estimation of conditional variance functions in stochastic regression*](https://doi.org/10.1093/biomet/85.3.645), Biometrika 85(3) (the residual-based estimator).
- J. Fan & J. Jiang (2000), [*Variable bandwidth and one-step local M-estimator*](https://doi.org/10.1007/BF02903849), Sci. China Ser. A 43 (the robust refinement).
"""
function conditional_variance(summary::Histogram2D, x::Real; h::Real,
                              kernel = gaussian_kernel, degree::Integer = 1, method::Symbol = :residual)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumY, sumY2 = _column_sums(summary)
  return _conditional_variance(cx, N, sumY, sumY2, x, h, kernel, degree, method)
end

function conditional_variance(summary::Histogram2D, xs::AbstractVector; h::Real,
                              kernel = gaussian_kernel, degree::Integer = 1, method::Symbol = :residual)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumY, sumY2 = _column_sums(summary)
  if method === :residual
    SR2 = _residual_response(cx, N, sumY, sumY2, h, kernel, degree)
    return [_clamp_var(local_fit(cx, N, SR2, x, h, kernel, degree)[1]) for x in xs]
  elseif method === :plugin
    return [_clamp_var(local_fit(cx, N, sumY2, x, h, kernel, degree)[1] -
                       local_fit(cx, N, sumY, x, h, kernel, degree)[1]^2) for x in xs]
  end
  throw(ArgumentError("conditional_variance `method` must be :residual (Fan–Yao, default) or " *
                      ":plugin (got :$method)"))
end

end # module ConditionalExpectations
