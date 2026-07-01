# ConditionalExpectations — §4.2: the conditional expectation E[g(Y)|x] of a smooth functional (the general
# reader), with the conditional mean m(x) = E[Y|x] and variance σ²(x) as its g = identity / y² cases, and
# the conditional pushforward density f_{g(Y)|X}. All by local-polynomial regression **in x** of a
# per-X-column **y-quadrature** response — so they compose the two atoms: kernel-weighting in x
# (`LocalPolynomial.local_fit`) ∘ power-sum quadrature in y (the `Expectations` coefficient primitive over
# the order-2 Y-power-sums). The estimand-first reference module — write every estimator docstring to this
# pattern (`docs/ARCHITECTURE.md` §2, §4.2).

module ConditionalExpectations

using LinearAlgebra: mul!

using ..Histograms: Histogram2D, _column_sums, _xbincenter, _ybincenter
using ..LocalPolynomial: local_fit, _leverage
using ..Kernels: gaussian_kernel
using ..Expectations: _interp_matrix

export conditional_mean, conditional_variance, conditional_expectation,
       conditional_expectation_standard_error, conditional_pushforward_density


# Local-polynomial order admitted by the conditional-moment estimators: 0 (Nadaraya–Watson),
# 1 (local-linear, the default), 2 (local-quadratic — exposes the curvature m'' that h⋆(x) needs).
@inline function _check_degree(degree)
  (degree == 0 || degree == 1 || degree == 2) ||
    throw(ArgumentError("degree must be 0, 1, or 2 (got degree = $degree)"))
end


# The per-X-column y-quadrature of a general smooth g: returns (cx, N, sumG), where for each X-column ix
# sumG[ix] = Σ over its Y-cells of the degree-2 local quadrature of g about the Y-bin center —
# a_0 n + a_1 S_01 + a_2 S_02, with (a_0, a_1, a_2) the coefficients of g's degree-2 interpolant on the
# Y-bin (order-2 is all `Histogram2D` stores in y, so this is exact for g of degree ≤ 2 — g = y gives
# Σ Y, g = y² gives Σ Y², matching `_column_sums` — and O(Δ_y²) for a general g). g's coefficients depend
# only on the Y-bin center, so they are formed once per Y-bin via the shared `Expectations._interp_matrix`.
function _column_response(summary::Histogram2D, g)
  counts = Array(summary.counts)
  S = Array(summary.powersums)
  nx = summary.nx
  ny = summary.ny
  halfwidth = 0.5 / summary.yscale                         # Δ_y/2
  nodes, vandermonde_inv = _interp_matrix(2)               # degree-2: the order-2 Y-power-sums
  gsamples = Vector{Float64}(undef, 3)
  bcoeffs = Vector{Float64}(undef, 3)
  a0 = Vector{Float64}(undef, ny)                          # g's degree-2 coefficients per Y-bin (in u = y − c^y)
  a1 = similar(a0)
  a2 = similar(a0)
  @inbounds for iy in 1:ny
    cy = _ybincenter(summary, iy)
    for m in 1:3
      gsamples[m] = g(cy + nodes[m] * halfwidth)
    end
    mul!(bcoeffs, vandermonde_inv, gsamples)               # coefficients in s = u/halfwidth
    a0[iy] = bcoeffs[1]
    a1[iy] = bcoeffs[2] / halfwidth                        # s → u
    a2[iy] = bcoeffs[3] / halfwidth^2
  end
  cx = [_xbincenter(summary, ix) for ix in 1:nx]
  N = zeros(Float64, nx)
  sumG = zeros(Float64, nx)
  @inbounds for iy in 1:ny
    c0 = a0[iy]; c1 = a1[iy]; c2 = a2[iy]
    for ix in 1:nx
      cell = ix + (iy - 1) * nx
      nc = counts[cell]
      iszero(nc) && continue
      N[ix] += nc
      sumG[ix] += c0 * nc + c1 * S[2, cell] + c2 * S[5, cell]   # a_0 n + a_1 S_01 + a_2 S_02
    end
  end
  return cx, N, sumG
end


"""
    conditional_expectation(summary::Histogram2D, g, x;
                            h, kernel = gaussian_kernel, degree = 1)   # ≡ ⟨g(Y) ∣ x⟩

**Estimand.** the conditional expectation of a smooth functional of the response,
``\\mathbb{E}[g(Y) \\mid X = x] = \\int g(y)\\, f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y``. For
``g = \\mathrm{identity}`` this is the regression function ``m(x)``
([`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean)); for
``g(y) = y^2`` it is ``\\mathbb{E}[Y^2 \\mid x]``, the plug-in ingredient of
[`conditional_variance`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_variance).

**Estimator.** a local-polynomial regression **in ``x``** of a per-``X``-column **``y``-quadrature** of
``g``. For each ``X``-column the response ``\\sum_i g(Y_i)`` is recovered from the stored order-2
``Y``-power-sums by the same quadrature the marginal
[`expectation`](@ref GBGStreamingStatistics.Expectations.expectation) uses
(``a_0 n + a_1 S_{01} + a_2 S_{02}``, with ``a_p`` the coefficients of ``g``'s degree-2 interpolant on the
``Y``-bin), then smoothed across ``x`` by [`local_fit`](@ref) of order `degree`; the estimate is the fitted
value ``\\beta_0``. It thus **composes the two atoms** — kernel-weighting in ``x`` after power-sum
quadrature in ``y`` — and reduces to
[`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean) at
``g = \\mathrm{identity}``.

Because `Histogram2D` stores only the order-2 ``Y``-power-sums, the ``y``-quadrature is **exact for ``g``
of degree ``\\le 2`` in ``y``** (so ``\\mathbb{E}[Y \\mid x]`` and ``\\mathbb{E}[Y^2 \\mid x]`` are exact)
and ``O(\\Delta_y^2)`` for a general smooth ``g`` — the price of the fixed order-2 joint summary. The
``x``-side carries the local-polynomial bias/variance of
[`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean). A vector of
query points returns a vector; where no data is local to ``x`` the estimate is `NaN`.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `g` — the response functional to average: any callable ``y \\mapsto g(y)`` (smooth on the ``Y``-bin
  scale). `g = identity` recovers the conditional mean, `g = abs2` the conditional ``\\mathbb{E}[Y^2 \\mid x]``.
- `x` — the conditioning value(s) at which the expectation is evaluated (scalar or vector).
- `h` — the conditioning (``x``) bandwidth (``h > 0``), supplied.
- `kernel` — the conditioning (``x``) kernel (any ``\\int K = 1`` callable; default `gaussian_kernel`).
- `degree` — the local-polynomial order in ``x``: `0` (Nadaraya–Watson), `1` (local-linear, the default),
  or `2`.
"""
function conditional_expectation(summary::Histogram2D, g, x::Real; h::Real,
                                 kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumG = _column_response(summary, g)
  return local_fit(cx, N, sumG, x, h, kernel, degree)[1]
end

function conditional_expectation(summary::Histogram2D, g, xs::AbstractVector; h::Real,
                                 kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumG = _column_response(summary, g)
  return [local_fit(cx, N, sumG, x, h, kernel, degree)[1] for x in xs]
end


# The standard error at one x from the pre-reduced column responses: √(leverage · Var(g(Y)|x)), with the
# plug-in conditional variance Ê[g²|x] − Ê[g|x]² and the local-regression leverage Σ_j N_j w_j². NaN where
# the fit is degenerate (no local data).
function _ce_standard_error(cx, N, sumG, sumG2, x, h, kernel, degree)
  eg = local_fit(cx, N, sumG, x, h, kernel, degree)[1]
  lev = _leverage(cx, N, x, h, kernel, degree)
  (isfinite(eg) && isfinite(lev)) || return oftype(eg, NaN)
  eg2 = local_fit(cx, N, sumG2, x, h, kernel, degree)[1]
  variance = eg2 - eg^2
  return sqrt(max(variance, zero(variance)) * max(lev, zero(lev)))
end


"""
    conditional_expectation_standard_error(summary::Histogram2D, g, x;
                                           h, kernel = gaussian_kernel, degree = 1)

**Estimand.** the sampling standard error of the estimator [`conditional_expectation`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation),
``\\mathrm{SE}\\bigl(\\widehat{\\mathbb{E}}[g(Y)\\mid x]\\bigr) =
\\sqrt{\\ell(x)\\,\\operatorname{Var}(g(Y)\\mid x)}``, with the local-regression **leverage**
``\\ell(x) = \\sum_j N_j\\, w_j(x)^2`` (the sum of squared effective weights of the fit's ``\\beta_0``)
and the conditional variance of the transformed response
``\\operatorname{Var}(g(Y)\\mid x) = \\mathbb{E}[g(Y)^2\\mid x] - \\mathbb{E}[g(Y)\\mid x]^2``.

**Estimator.** the **leverage × noise-scale** readout
``\\sqrt{\\hat\\ell(x)\\,(\\widehat{\\mathbb{E}}[g^2\\mid x] - \\widehat{\\mathbb{E}}[g\\mid x]^2)}``: the
leverage ``\\hat\\ell(x)`` from the local-fit design (its plain and **squared-kernel** moments; read off the
actual weights, so correct at any `degree`) times the plug-in conditional variance of ``g(Y)`` — itself two
[`conditional_expectation`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_expectation)
readouts (``\\mathbb{E}[g^2\\mid x]`` the ``g^2`` reader). This is the exact analog of the marginal
[`expectation_standard_error`](@ref GBGStreamingStatistics.Expectations.expectation_standard_error)
``= \\sqrt{\\operatorname{Var}(g(X))\\,/\\,n}``, with the global-mean leverage ``1/n`` replaced by the local
``\\sum_j N_j w_j^2``; at ``g = \\mathrm{identity}`` it is the standard error of the conditional mean
``m(x)``. It is the plug-in sampling error (not the ``O(h^2)`` smoothing bias). `NaN` where no data is
local to ``x``.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `g` — the response functional whose conditional average's standard error is wanted (as in
  [`conditional_expectation`](@ref
  GBGStreamingStatistics.ConditionalExpectations.conditional_expectation)); `g = identity` gives the
  conditional mean's standard error.
- `x` — the conditioning value(s) (scalar or vector).
- `h` — the conditioning (``x``) bandwidth (``h > 0``), supplied.
- `kernel` — the conditioning (``x``) kernel (any ``\\int K = 1`` callable; default `gaussian_kernel`).
- `degree` — the local-polynomial order in ``x``: `0`, `1` (the default), or `2`.
"""
function conditional_expectation_standard_error(summary::Histogram2D, g, x::Real; h::Real,
                                                kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumG = _column_response(summary, g)
  sumG2 = _column_response(summary, y -> g(y)^2)[3]
  return _ce_standard_error(cx, N, sumG, sumG2, x, h, kernel, degree)
end

function conditional_expectation_standard_error(summary::Histogram2D, g, xs::AbstractVector; h::Real,
                                                kernel = gaussian_kernel, degree::Integer = 1)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  cx, N, sumG = _column_response(summary, g)
  sumG2 = _column_response(summary, y -> g(y)^2)[3]
  return [_ce_standard_error(cx, N, sumG, sumG2, x, h, kernel, degree) for x in xs]
end


# conditional_mean is the g = identity case of conditional_expectation, folded in exactly (the y-quadrature
# of identity is Σ Y, degree 1 ≤ 2 so exact); x may be a scalar or a vector (dispatched by the general one).
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
conditional_mean(summary::Histogram2D, x; h::Real, kernel = gaussian_kernel, degree::Integer = 1) =
  conditional_expectation(summary, identity, x; h = h, kernel = kernel, degree = degree)


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


"""
    conditional_pushforward_density(summary::Histogram2D, g, x, z;
                                    hx, hz, kernelx = gaussian_kernel, kernelz = gaussian_kernel,
                                    degree = 1)   # ≡ f_{g(Y) ∣ X}

**Estimand.** the density of the transformed response ``Z = g(Y)`` given ``X = x`` — the **conditional
pushforward**, change of variables applied to the conditional law,
``f_{g(Y)\\mid X}(z \\mid x) = \\sum_{y:\\,g(y)=z} f_{Y\\mid X}(y \\mid x)\\,/\\,|g'(y)|``.

**Estimator.** the conditional kernel density of ``g(Y)``, written as a conditional expectation and
evaluated by [`conditional_expectation`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation):

```math
\\hat f_{g(Y)\\mid X}(z \\mid x) \\;=\\; \\widehat{\\mathbb{E}}\\bigl[\\, K_{h_z}(g(Y) - z) \\mid x \\,\\bigr] ,
```

i.e. `conditional_expectation(summary, y -> K_hz(g(y) - z), x)`. No inverse ``g^{-1}`` and no explicit
Jacobian is formed; a **non-monotone** ``g`` contributes from every preimage automatically, ``g'=0`` is an
integrable caustic, and a symmetric response kernel ``K_z`` (`kernelz`, ``\\int K = 1``) makes
``\\int \\hat f\\,\\mathrm dz = 1`` exactly. Two smoothing directions: the ``x``-conditioning regression
(`hx`, `kernelx`, `degree` — **degree-aware** in ``x``, unlike the degree-0 double-kernel
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)) and the ``z``
smoothing (`hz`, `kernelz`). At ``g = \\mathrm{identity}`` with `degree = 0` it coincides with
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density) — which stays
the specialized native estimator (with the bounded-response asymmetric `kernely` and `normalize`); use it
directly for the untransformed response. Needs ``g`` smooth on the ``Y``-bin scale and the response
bandwidth ``h_z \\gtrsim |g'|\\,\\Delta_y``. A vector `z` returns the density curve at the fixed
conditioning value `x`.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `g` — the transform defining ``Z = g(Y)`` (any callable, smooth on the ``Y``-bin scale).
- `x` — the conditioning value (a scalar).
- `z` — the transformed-response point(s) at which the density is evaluated (a scalar, or a vector for the
  curve).
- `hx`, `hz` — the conditioning (``x``) and transformed-response (``z``) bandwidths (both ``> 0``), supplied.
- `kernelx` — the conditioning (``x``) kernel (any ``\\int K = 1`` callable; default `gaussian_kernel`).
- `kernelz` — the response-direction (``z``) kernel (default `gaussian_kernel`); a symmetric kernel
  integrates to 1 exactly.
- `degree` — the local-polynomial order in ``x``: `0`, `1` (the default), or `2`.
"""
function conditional_pushforward_density(summary::Histogram2D, g, x::Real, z::Real;
                                         hx::Real, hz::Real, kernelx = gaussian_kernel,
                                         kernelz = gaussian_kernel, degree::Integer = 1)
  hz > 0 || throw(ArgumentError("response bandwidth hz must be positive (got hz = $hz)"))
  return conditional_expectation(summary, y -> kernelz((g(y) - z) / hz) / hz, x;
                                 h = hx, kernel = kernelx, degree = degree)
end

function conditional_pushforward_density(summary::Histogram2D, g, x::Real, zs::AbstractVector;
                                         hx::Real, hz::Real, kernelx = gaussian_kernel,
                                         kernelz = gaussian_kernel, degree::Integer = 1)
  hz > 0 || throw(ArgumentError("response bandwidth hz must be positive (got hz = $hz)"))
  return [conditional_expectation(summary, y -> kernelz((g(y) - z) / hz) / hz, x;
                                  h = hx, kernel = kernelx, degree = degree) for z in zs]
end

end # module ConditionalExpectations
