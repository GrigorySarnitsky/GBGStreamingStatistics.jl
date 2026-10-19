# ConditionalExpectations — §4.2: the conditional expectation 𝔼[g(Y) | X = x] of a smooth functional (the general
# reader), with the conditional mean m(x) = 𝔼[Y | X = x] and variance σ²(x) as its g = identity / y² cases, and
# the conditional pushforward density f_{g(Y)|X}. All by local-polynomial regression in x of a
# per-x-column y-quadrature response — so they compose the two atoms: kernel-weighting in x
# (`LocalPolynomial.local_fit`) ∘ power-sum quadrature in y (the `Expectations` coefficient primitive over
# the order-2 y-power-sums). The estimand-first reference module — write every estimator docstring to this
# pattern (`docs/ARCHITECTURE.md` §2, §4.2).

module ConditionalExpectations

using LinearAlgebra: mul!

using ..Histograms: Histogram2D, _column_sums, _xbincenter, _ybincenter, _total
using ..LocalPolynomial: local_fit, _leverage
using ..Kernels: gaussian_kernel, _quad, AsymmetricKernel
using ..Expectations: _interp_matrix

export conditional_mean, conditional_variance, conditional_variance_standard_error,
       conditional_variance_with_standard_error, conditional_expectation,
       conditional_expectation_standard_error, conditional_expectation_with_standard_error,
       conditional_pushforward_density, conditional_pushforward_density_standard_error,
       conditional_pushforward_density_with_standard_error


# Local-polynomial order admitted by the conditional-moment estimators: 0 (Nadaraya–Watson),
# 1 (local-linear, the default), 2 (local-quadratic — exposes the curvature m̂'' that h⋆(x) needs).
@inline function _check_degree(degree)
  (degree == 0 || degree == 1 || degree == 2) ||
    throw(ArgumentError("degree must be 0, 1, or 2 (got degree = $degree)"))
end


# The per-x-column y-quadrature of a general smooth g: returns (c_x, N, sum_G), where for each x-column i_x
# sum_G[i_x] = Σ over its y-cells of the degree-2 local quadrature of g about the y-bin center —
# a_0 n_{jl} + a_1 S_01 + a_2 S_02, with (a_0, a_1, a_2) the coefficients of g's degree-2 interpolant on the
# y-bin (order-2 is all `Histogram2D` stores in y, so this is exact for g of degree ≤ 2 — g = y gives
# Σ Y, g = y² gives Σ Y², matching `_column_sums` — and O(Δ_y²) for a general g). g's coefficients depend
# only on the y-bin center, so they are formed once per y-bin via the shared `Expectations._interp_matrix`.
function _column_response(summary::Histogram2D, g)
  counts = Array(summary.counts)
  S = Array(summary.powersums)
  n_x = summary.n_x
  n_y = summary.n_y
  halfwidth = 0.5 / summary.scale_y                         # Δ_y/2
  nodes, vandermonde_inv = _interp_matrix(2)               # degree-2: the order-2 y-power-sums
  gsamples = Vector{Float64}(undef, 3)
  bcoeffs = Vector{Float64}(undef, 3)
  a0 = Vector{Float64}(undef, n_y)                          # g's degree-2 coefficients per y-bin (in y − c^y)
  a1 = similar(a0)
  a2 = similar(a0)
  @inbounds for i_y in 1:n_y
    c_y = _ybincenter(summary, i_y)
    for m in 1:3
      gsamples[m] = g(c_y + nodes[m] * halfwidth)
    end
    mul!(bcoeffs, vandermonde_inv, gsamples)               # coefficients in s = (y − c^y)/halfwidth
    a0[i_y] = bcoeffs[1]
    a1[i_y] = bcoeffs[2] / halfwidth                        # s → y − c^y
    a2[i_y] = bcoeffs[3] / halfwidth^2
  end
  c_x = [_xbincenter(summary, i_x) for i_x in 1:n_x]
  N = zeros(Float64, n_x)
  sum_G = zeros(Float64, n_x)
  @inbounds for i_y in 1:n_y
    c0 = a0[i_y]; c1 = a1[i_y]; c2 = a2[i_y]
    for i_x in 1:n_x
      cell = i_x + (i_y - 1) * n_x
      nc = counts[cell]
      iszero(nc) && continue
      N[i_x] += nc
      sum_G[i_x] += c0 * nc + c1 * S[2, cell] + c2 * S[5, cell]   # a_0 n_{jl} + a_1 S_01 + a_2 S_02
    end
  end
  return c_x, N, sum_G
end


"""
    conditional_expectation(summary::Histogram2D, g, x;
                            kernel = gaussian_kernel,
                            degree = kernel isa AsymmetricKernel ? 0 : 1, h)

**Estimand.** The estimand is the conditional expectation of a smooth functional of the response,
``\\mathbb{E}[g(Y) \\mid X = x] = \\int g(y)  f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y``. For
``g = \\mathrm{identity}`` this is the regression function ``m(x)``
([`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean)); for
``g(y) = y^2`` it is ``\\mathbb{E}[Y^2 \\mid X = x]``, the plug-in ingredient of
[`conditional_variance`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_variance).

**Estimator.** The estimate is a local-polynomial regression in ``x`` of a per-``x``-column (per conditioning bin)
``y``-quadrature of ``g``. Per ``x``-column the response ``\\sum_i g(Y_i)`` is
recovered from the stored order-2 ``y``-power-sums (per-cell sums of the response and its square; exact
for ``g`` of degree ``\\le 2``) by the marginal
[`expectation`](@ref GBGStreamingStatistics.Expectations.expectation)'s quadrature, then smoothed across ``x`` by a
local-polynomial fit of order `degree`; the estimate is the fitted value ``\\hat\\beta_0``. It thus applies
kernel-weighting in ``x`` after power-sum quadrature in ``y``, and reduces to `conditional_mean` at
``g = \\mathrm{identity}``.

Because a [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) stores only the order-2
``y``-power-sums, the ``y``-quadrature is exact for ``g`` of degree ``\\le 2`` in ``y`` (so
``\\mathbb{E}[Y \\mid X = x]`` and ``\\mathbb{E}[Y^2 \\mid X = x]`` are exact) and ``O(\\Delta_y^2)`` in the response-bin width for a general
smooth ``g``. The ``x``-side carries the local-polynomial
bias/variance of `conditional_mean`. A vector of query points returns a vector, `NaN` where no data is local to
``x``.

# Arguments

- `summary` — the joint `Histogram2D` of ``(X, Y)``.
- `g` — the response functional to average: any callable ``y \\mapsto g(y)`` (well approximated by a
  quadratic within each ``y``-bin). `g = identity` gives the conditional mean, `g = abs2` the conditional ``\\mathbb{E}[Y^2 \\mid X = x]``.
- `x` — the conditioning value(s) where the expectation is evaluated (scalar or vector).
- `kernel` — the conditioning (``x``) kernel: any callable kernel of unit integral, ``\\int K = 1`` (default
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)) or an
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel)
  ([`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) /
  [`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel), for bounded conditioning).
- `degree` — the local-polynomial order in ``x``: `0` (Nadaraya–Watson), `1` (local-linear), or `2`. The
  default is `1`, or `0` under an asymmetric `kernel` (which fits no polynomial).
- `h` — the conditioning (``x``) bandwidth (``h > 0``), supplied.

Reference: J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its Applications*](https://doi.org/10.1201/9780203748725),
Chapman & Hall.
"""
function conditional_expectation(summary::Histogram2D, g, x::Real;
                                 kernel = gaussian_kernel,
                                 degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                 h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  return local_fit(c_x, N, sum_G, x, kernel, degree, h)[1]
end

function conditional_expectation(summary::Histogram2D, g, xs::AbstractVector;
                                 kernel = gaussian_kernel,
                                 degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                 h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  return [local_fit(c_x, N, sum_G, x, kernel, degree, h)[1] for x in xs]
end


# The value and its standard error at one x from the pre-reduced column responses: the local-fit β̂₀ and
# √(leverage · Var(g(Y) | X = x)), with the plug-in conditional variance 𝔼̂[g² | X = x] − 𝔼̂[g | X = x]² and
# the leverage Σ_j N_j ω_j². Both NaN where the fit is degenerate (no local data). Shared by the standard-error reader
# and the paired value-plus-error reader, so the two agree by construction.
function _ce_value_and_se(c_x, N, sum_G, sum_G2, x, kernel, degree, h)
  eg = local_fit(c_x, N, sum_G, x, kernel, degree, h)[1]
  lev = _leverage(c_x, N, x, kernel, degree, h)
  (isfinite(eg) && isfinite(lev)) || return (eg, oftype(eg, NaN))
  eg2 = local_fit(c_x, N, sum_G2, x, kernel, degree, h)[1]
  variance = eg2 - eg^2
  return (eg, sqrt(max(variance, zero(variance)) * max(lev, zero(lev))))
end


"""
    conditional_expectation_standard_error(summary::Histogram2D, g, x;
                                           kernel = gaussian_kernel,
                                           degree = kernel isa AsymmetricKernel ? 0 : 1,
                                           h)

**Estimand.** The estimand is the sampling standard error of the estimator `conditional_expectation`,
``\\operatorname{SE}\\bigl(\\widehat{\\mathbb{E}}[g(Y)\\mid X = x]\\bigr) =
\\sqrt{\\ell(x) \\operatorname{Var}(g(Y)\\mid X = x)}``, with the local-regression *leverage*
``\\ell(x) = \\sum_j N_j  \\omega_j(x)^2`` (``N_j`` the per-``x``-column count, ``\\omega_j`` the effective weights of the fit's ``\\hat\\beta_0``)
and the conditional variance of the transformed response
``\\operatorname{Var}(g(Y)\\mid X = x) = \\mathbb{E}[g(Y)^2\\mid X = x] - \\mathbb{E}[g(Y)\\mid X = x]^2``.

**Estimator.** The estimate is the leverage times the noise scale,
``\\sqrt{\\hat\\ell(x) (\\widehat{\\mathbb{E}}[g^2\\mid X = x] - \\widehat{\\mathbb{E}}[g\\mid X = x]^2)}``: the
leverage ``\\hat\\ell(x)`` of the local fit (correct at any `degree`) times the plug-in conditional variance
of ``g(Y)`` — itself two evaluations of `conditional_expectation` (``\\mathbb{E}[g^2\\mid X = x]`` with ``g^2`` for the
functional). This is the analog of the marginal
[`expectation_standard_error`](@ref GBGStreamingStatistics.Expectations.expectation_standard_error),
``\\sqrt{\\operatorname{Var}(g(X)) / n}``, with the global-mean leverage ``1/n`` replaced by the local
``\\sum_j N_j \\omega_j^2``; at ``g = \\mathrm{identity}`` it is the standard error of the conditional mean
``m(x)``. It is the plug-in sampling error (not the ``O(h^2)`` smoothing bias). It is `NaN` where no data
is local to ``x``.

Arguments: as for [`conditional_expectation`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_expectation).

"""
function conditional_expectation_standard_error(summary::Histogram2D, g, x::Real;
                                                kernel = gaussian_kernel,
                                                degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  sum_G2 = _column_response(summary, y -> g(y)^2)[3]
  return _ce_value_and_se(c_x, N, sum_G, sum_G2, x, kernel, degree, h)[2]
end

function conditional_expectation_standard_error(summary::Histogram2D, g, xs::AbstractVector;
                                                kernel = gaussian_kernel,
                                                degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  sum_G2 = _column_response(summary, y -> g(y)^2)[3]
  return [_ce_value_and_se(c_x, N, sum_G, sum_G2, x, kernel, degree, h)[2] for x in xs]
end


"""
    conditional_expectation_with_standard_error(summary::Histogram2D, g, x;
        kernel = gaussian_kernel, degree = kernel isa AsymmetricKernel ? 0 : 1, h)

Returns the estimate and its standard error from two column reductions (the functional, then its
square) on one local-polynomial design — the `NamedTuple`
``(\\texttt{value}, \\texttt{standard\\_error})``, with `value` the `conditional_expectation`
``\\widehat{\\mathbb{E}}[g(Y)\\mid X = x]`` and `standard_error` its
[`conditional_expectation_standard_error`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation_standard_error). They describe the
same estimator by construction, their ``y``-quadrature responses fitted on one design. For a vector `x`,
returns a `NamedTuple` of vectors.

Arguments: as for [`conditional_expectation`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation).

"""
function conditional_expectation_with_standard_error(summary::Histogram2D, g, x::Real;
                                                     kernel = gaussian_kernel,
                                                     degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                     h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  sum_G2 = _column_response(summary, y -> g(y)^2)[3]
  v, se = _ce_value_and_se(c_x, N, sum_G, sum_G2, x, kernel, degree, h)
  return (value = v, standard_error = se)
end

function conditional_expectation_with_standard_error(summary::Histogram2D, g, xs::AbstractVector;
                                                     kernel = gaussian_kernel,
                                                     degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                     h::Real)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_G = _column_response(summary, g)
  sum_G2 = _column_response(summary, y -> g(y)^2)[3]
  vses = [_ce_value_and_se(c_x, N, sum_G, sum_G2, x, kernel, degree, h) for x in xs]
  return (value = [vs[1] for vs in vses], standard_error = [vs[2] for vs in vses])
end


# conditional_mean is the g = identity case of conditional_expectation, folded in exactly (the y-quadrature
# of identity is Σ Y, degree 1 ≤ 2 so exact); x may be a scalar or a vector (dispatched by the general one).
"""
    conditional_mean(summary::Histogram2D, x; kernel = gaussian_kernel,
                     degree = kernel isa AsymmetricKernel ? 0 : 1, h)

**Estimand.** The estimand is the conditional expectation
``m(x) = \\mathbb{E}[Y \\mid X = x] = \\int y  f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y`` — the
``L^2(\\mathbb{P}_X)`` projection of the response ``Y`` onto the conditioning variable ``X``, i.e. the
regression function.

**Estimator** (local-polynomial; Fan & Gijbels 1996). From the joint `summary` (a
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X, Y)``) form the design
and response sums, weighted by the kernel scaled to the bandwidth, ``U_p(x) = \\sum_i K_h(X_i-x)(X_i-x)^p`` and
``V_p(x) = \\sum_i K_h(X_i-x)(X_i-x)^p Y_i`` — each sample taken at its conditioning-bin center (an error first
order in the bin width), from the per-bin counts and exact response totals ``\\sum Y`` — and solve the local
normal equations of order `degree` for the *jet* ``(\\hat\\beta_0, \\dots, \\hat\\beta_{\\texttt{degree}})`` — the
local-polynomial coefficients at ``x`` (``\\hat\\beta_0`` the fitted value, ``\\hat\\beta_1`` the slope, ``2\\hat\\beta_2``
the curvature); the estimate is the fitted value ``\\hat m(x) = \\hat\\beta_0``. At the default `degree = 1`
(local-linear) this is

```math
\\hat m(x) = \\frac{U_2 V_0 - U_1 V_1}{U_0 U_2 - U_1^{2}},
```

which removes the design-density bias ``2 m'(x) f_X'(x)/f_X(x)`` (the design density being the marginal
density of the conditioning variable) of the Nadaraya–Watson reduction
``\\hat m = V_0/U_0`` (the `degree = 0` local-constant fit). The bandwidth ``h`` is supplied (or the
response-aware local plug-in ``h^\\star(x)`` from
[`bandwidth(summary, x; rule = :local)`](@ref GBGStreamingStatistics.Bandwidth.bandwidth));
`degree = 2` (local-quadratic) additionally exposes the curvature ``\\hat m''(x)`` that ``h^\\star(x)`` needs. A vector of query points
returns a vector. It is consistent, with leading-order ``\\operatorname{bias}\\hat m(x) \\approx
\\tfrac{1}{2} h^2 \\mu_2(K)  m''(x)`` and ``\\operatorname{Var}\\hat m(x) \\approx
R(K) \\sigma^2(x)/(n h f_X(x))``, in the kernel's second moment and roughness, the conditional variance
and the sample count.

**Boundary adaptivity.** The local-linear fit (`degree = 1`, the default) is *boundary-adaptive* (Fan &
Gijbels 1996): conditioning on a bounded or positive quantity — a droplet radius ``X \\ge 0``, say — it
is first-order correct near the edge with the symmetric `kernel`, correcting the ``O(h)`` boundary bias to
the interior ``O(h^2)`` automatically. The `degree = 0` Nadaraya–Watson fit keeps the ``O(h)`` edge bias,
so near a bounded edge use `degree = 1`, which is also unbiased for a linear mean. A `degree = 0`
asymmetric `kernel` ([`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) /
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel)) is *admitted* too, but these boundary
kernels are a density tool whose asymmetry interacts with the regression design density, so they do not
generally improve the conditional mean. A bounded *response* is the other case: there the response kernel
can be asymmetric —
see [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density).

# Arguments

- `summary` — the joint `Histogram2D` of ``(X, Y)``.
- `x` — the conditioning value(s) at which ``m`` is evaluated (a scalar for one point, or a vector
  for the whole curve).
- `kernel` — the smoothing kernel: a symmetric callable ``K`` (``\\int K = 1``; default
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel), at any `degree`), or an
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel) (`gamma_kernel` / `beta_kernel`) for a bounded/positive conditioning
  variable (`degree = 0` only).
- `degree` — the local-polynomial order: `0` (Nadaraya–Watson), `1` (local-linear), or `2`
  (local-quadratic). The default is `1`, or `0` under an asymmetric `kernel`, which fits no polynomial and admits only
  `0`.
- `h` — the smoothing bandwidth (``h > 0``), supplied.

Reference: J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its Applications*](https://doi.org/10.1201/9780203748725),
Chapman & Hall.
"""
conditional_mean(summary::Histogram2D, x; kernel = gaussian_kernel,
                 degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1), h::Real) =
  conditional_expectation(summary, identity, x; h = h, kernel = kernel, degree = degree)


# Clamp a variance estimate to be non-negative; preserve NaN (off-support / degenerate design).
@inline _clamp_var(v) = v > 0 ? v : oftype(v, isnan(v) ? NaN : 0.0)


# Fan & Yao (1998) residual response for the bias-reduced conditional variance: per x-column, the total
# of squared residuals Σ(Y − m̂(X))², with m̂(X) approximated by the column-center fit m̂(c^x_j)
# (exact to O(Δ), the bin width). The column total is then ΣY² − 2 m̂(c^x_j) ΣY + m̂(c^x_j)² N_j from the
# stored column sums. A second local fit of this response (≡ 𝔼̂[(Y − m̂)² | X = x]) estimates σ²(x) with the
# curved-𝔼[Y² | X = x] bias of the plug-in removed — 𝔼[(Y − m)² | X = x] = σ²(x) is far flatter than
# 𝔼[Y² | X = x].
function _residual_response(c_x, N, sum_Y, sum_Y2, kernel, degree, h)
  SR2 = similar(sum_Y2)
  @inbounds for j in eachindex(c_x)
    if N[j] > 0
      μj = local_fit(c_x, N, sum_Y, c_x[j], kernel, degree, h)[1]
      SR2[j] = sum_Y2[j] - 2 * μj * sum_Y[j] + μj^2 * N[j]
    else
      SR2[j] = zero(eltype(SR2))
    end
  end
  return SR2
end


# σ²(x) at a single point by the chosen method, from the pre-reduced column sums.
function _conditional_variance(c_x, N, sum_Y, sum_Y2, x, kernel, degree, h, method)
  if method === :residual
    SR2 = _residual_response(c_x, N, sum_Y, sum_Y2, kernel, degree, h)
    return _clamp_var(local_fit(c_x, N, SR2, x, kernel, degree, h)[1])
  elseif method === :plugin
    m = local_fit(c_x, N, sum_Y, x, kernel, degree, h)[1]
    ey2 = local_fit(c_x, N, sum_Y2, x, kernel, degree, h)[1]
    return _clamp_var(ey2 - m^2)
  end
  throw(ArgumentError("conditional_variance `method` must be :residual (Fan–Yao, default) or " *
                      ":plugin (got :$method)"))
end


"""
    conditional_variance(summary::Histogram2D, x; kernel = gaussian_kernel,
                         degree = kernel isa AsymmetricKernel ? 0 : 1, h,
                         method = :residual)

**Estimand.** The estimand is the conditional variance
``\\sigma^2(x) = \\operatorname{Var}(Y \\mid X = x) = \\mathbb{E}[Y^2 \\mid X = x] - m(x)^2``, with
``m(x)`` the [`conditional_mean`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_mean).

**Estimator** (`method = :residual`, the default — Fan & Yao 1998). The two-stage residual-based
estimator: fit ``\\hat m`` at each x-column center, then take a second local fit of the per-column
squared residuals, ``\\hat\\sigma^2(x) = \\widehat{\\mathbb{E}}\\bigl[(Y - \\hat m(X))^2 \\mid X = x\\bigr]``.
Because ``\\mathbb{E}[(Y - m(X))^2 \\mid X = x] = \\sigma^2(x)`` is far flatter than
``\\mathbb{E}[Y^2 \\mid X = x]``, its local fit is much less biased — it removes the local-polynomial
smoothing bias of the curved ``\\mathbb{E}[Y^2 \\mid X = x]`` that the plug-in carries. The
squared-residual total of ``x``-column ``j`` is recovered exactly from its stored response sums and sample count,
``\\sum Y^2 - 2\\hat m \\sum Y + \\hat m^2 N_j``.

**Estimator** (`method = :plugin`). The estimate is the difference ``\\widehat{\\mathbb{E}}[Y^2 \\mid X = x] -
\\hat m(x)^2`` (the same local fit applied to the ``Y^2`` response), clamped non-negative; cheaper
(no first stage), it carries the curved-``\\mathbb{E}[Y^2 \\mid X = x]`` smoothing bias. A vector `x` returns a
vector, `NaN` where no data is local.

# Arguments

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``\\sigma^2`` is evaluated (scalar or vector).
- `kernel` — the smoothing kernel ``K`` (any callable with ``\\int K = 1``; default
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)), or an
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel) (`degree = 0` only).
- `degree` — the local-polynomial order for the underlying fits: `0`, `1`, or `2`. The default is `1`, or
  `0` under an asymmetric `kernel`, which fits no polynomial and admits only `0`.
- `h` — the smoothing bandwidth (``h > 0``), supplied.
- `method` — `:residual` (the default, the bias-reduced Fan & Yao two-stage estimator) or `:plugin`
  (the simple ``\\widehat{\\mathbb{E}}[Y^2 \\mid X = x] - \\hat m^2`` difference).

# References

- J. Fan & Q. Yao (1998), [*Efficient estimation of conditional variance functions in stochastic regression*](https://doi.org/10.1093/biomet/85.3.645), Biometrika 85(3) (the residual-based estimator).
- J. Fan & J. Jiang (2000), [*Variable bandwidth and one-step local M-estimator*](https://doi.org/10.1007/BF02903849), Sci. China Ser. A 43 (a robust one-step M-estimator refinement, not offered).
"""
function conditional_variance(summary::Histogram2D, x::Real; kernel = gaussian_kernel,
                              degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                              h::Real, method::Symbol = :residual)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_Y, sum_Y2 = _column_sums(summary)
  return _conditional_variance(c_x, N, sum_Y, sum_Y2, x, kernel, degree, h, method)
end

function conditional_variance(summary::Histogram2D, xs::AbstractVector;
                              kernel = gaussian_kernel,
                              degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                              h::Real, method::Symbol = :residual)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  _check_degree(degree)
  c_x, N, sum_Y, sum_Y2 = _column_sums(summary)
  if method === :residual
    SR2 = _residual_response(c_x, N, sum_Y, sum_Y2, kernel, degree, h)
    return [_clamp_var(local_fit(c_x, N, SR2, x, kernel, degree, h)[1]) for x in xs]
  elseif method === :plugin
    return [_clamp_var(local_fit(c_x, N, sum_Y2, x, kernel, degree, h)[1] -
                       local_fit(c_x, N, sum_Y, x, kernel, degree, h)[1]^2) for x in xs]
  end
  throw(ArgumentError("conditional_variance `method` must be :residual (Fan–Yao, default) or " *
                      ":plugin (got :$method)"))
end


"""
    conditional_variance_standard_error(summary::Histogram2D, x;
                                        kernel = gaussian_kernel,
                                        degree = kernel isa AsymmetricKernel ? 0 : 1,
                                        h)

**Estimand.** The estimand is the sampling standard error of the conditional variance
`conditional_variance`,
``\\operatorname{SE}(\\hat\\sigma^2(x)) = \\sqrt{\\ell(x) \\bigl(\\mu_4(x) - \\sigma^4(x)\\bigr)}`` with the
local-regression leverage ``\\ell(x)``, the conditional variance ``\\sigma^2(x)``, and the 4th central
conditional moment ``\\mu_4(x) = \\mathbb{E}[(Y - m(x))^4 \\mid X = x]``.

**Estimator.** The residual approximation: since ``\\hat\\sigma^2(x)`` is a conditional mean of the
squared residual ``(Y - \\hat m(x))^2``, its standard error is the
[`conditional_expectation_standard_error`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation_standard_error) of that response,
``\\sqrt{\\hat\\ell(x) (\\hat\\mu_4(x) - \\hat\\sigma^4(x))}``, with ``\\hat\\mu_4`` the ``y``-quadrature of the
smooth ``(Y - \\hat m(x))^4`` off the order-2 ``y``-power-sums — so the quadrature, ``O(\\Delta_y^2)`` in the response-bin width,
captures the conditional kurtosis (a heavy-tailed response widens the band) with no 4th-moment storage.
It is the leading-order SE common to both variance methods. A vector `x` returns a vector, `NaN` where
no data is local.

Arguments: as for [`conditional_variance`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_variance), without `method` (the leading-order
SE is common to `:residual` and `:plugin`).

"""
function conditional_variance_standard_error(summary::Histogram2D, x::Real;
                                             kernel = gaussian_kernel,
                                             degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                             h::Real)
  m = conditional_mean(summary, x; h = h, kernel = kernel, degree = degree)
  return conditional_expectation_standard_error(summary, y -> (y - m)^2, x; h = h, kernel = kernel,
                                                degree = degree)
end

function conditional_variance_standard_error(summary::Histogram2D, xs::AbstractVector;
                                             kernel = gaussian_kernel,
                                             degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                             h::Real)
  return [conditional_variance_standard_error(summary, x; h = h, kernel = kernel, degree = degree) for x in xs]
end


"""
    conditional_variance_with_standard_error(summary::Histogram2D, x;
                                             kernel = gaussian_kernel,
                                             degree = kernel isa AsymmetricKernel ? 0 : 1,
                                             h, method = :residual)

Returns the conditional variance
`conditional_variance` and its
standard error as the `NamedTuple` ``(\\texttt{value}, \\texttt{standard\\_error})``, with `value` the
variance (`method = :residual`, the default, or `:plugin`) and `standard_error` its
[`conditional_variance_standard_error`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_variance_standard_error) (the leading-order
residual approximation, common to both methods). For a vector `x` it returns a `NamedTuple` of vectors.

Arguments: as for [`conditional_variance`](@ref GBGStreamingStatistics.ConditionalExpectations.conditional_variance).

"""
function conditional_variance_with_standard_error(summary::Histogram2D, x::Real;
                                                  kernel = gaussian_kernel,
                                                  degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                  h::Real, method::Symbol = :residual)
  v = conditional_variance(summary, x; h = h, kernel = kernel, degree = degree, method = method)
  se = conditional_variance_standard_error(summary, x; h = h, kernel = kernel, degree = degree)
  return (value = v, standard_error = se)
end

function conditional_variance_with_standard_error(summary::Histogram2D, xs::AbstractVector;
                                                  kernel = gaussian_kernel,
                                                  degree::Integer = (kernel isa AsymmetricKernel ? 0 : 1),
                                                  h::Real, method::Symbol = :residual)
  return (value = conditional_variance(summary, xs; h = h, kernel = kernel, degree = degree, method = method),
          standard_error = conditional_variance_standard_error(summary, xs; h = h, kernel = kernel, degree = degree))
end


"""
    conditional_pushforward_density(summary::Histogram2D, g, x, z;
                                    kernel_x = gaussian_kernel,
                                    kernel_z = gaussian_kernel,
                                    degree = kernel_x isa AsymmetricKernel ? 0 : 1, h_x,
                                    h_z)

**Estimand.** The estimand is the density of the response transformed by a smooth function, ``Z = g(Y)``, given ``X = x`` — the *conditional
pushforward*, change of variables applied to the conditional law,
``f_{g(Y)\\mid X}(z \\mid x) = \\sum_{y: g(y)=z} f_{Y\\mid X}(y \\mid x) / |g'(y)|``.

**Estimator.** The estimate is the conditional kernel density of ``g(Y)``, written as a conditional
expectation and evaluated by [`conditional_expectation`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_expectation):

```math
\\hat f_{g(Y)\\mid X}(z \\mid x) = \\widehat{\\mathbb{E}}\\bigl[  K_{h_z}(g(Y) - z) \\mid X = x  \\bigr] ,
```

with the response kernel scaled to the bandwidth `h_z`. No inverse ``g^{-1}`` or Jacobian is formed: a
non-monotone ``g`` contributes from every preimage, ``g'=0`` (where the Jacobian diverges) is integrable,
and a symmetric response kernel ``K_z`` (`kernel_z`, ``\\int K = 1``) makes
``\\int \\hat f_{g(Y)\\mid X}\\,\\mathrm dz = 1`` exactly. It smooths in ``x`` (`h_x`, `kernel_x`, `degree` —
degree-aware in ``x``, unlike the degree-0
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)) and in ``z``
(`h_z`, `kernel_z`). At ``g = \\mathrm{identity}`` with `degree = 0` it matches `conditional_density`
(which ignores the within-cell spread of the response) to second order in the response-bin width; the
latter, with its asymmetric response kernel and exact normalization, stays native for the untransformed
response. It needs ``g`` near-quadratic within each ``y``-bin and ``h_z \\gtrsim |g'| \\Delta_y``, the bin
width stretched by the slope.

# Arguments

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `g` — the transform defining ``Z = g(Y)`` (any callable, near-quadratic within each ``y``-bin).
- `x` — the conditioning value (a scalar).
- `z` — the transformed-response point(s) at which the density is evaluated (scalar, or a vector for the
  curve).
- `kernel_x` — the conditioning (``x``) kernel (a callable kernel with ``\\int K = 1``, or an
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel); default
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)).
- `kernel_z` — the response-direction (``z``) kernel (default `gaussian_kernel`; a symmetric one integrates
  to 1 exactly).
- `degree` — the local-polynomial order in ``x``: `0`, `1`, or `2`. The default is `1`, or `0` under an
  asymmetric `kernel_x`, which fits no polynomial and admits only `0`.
- `h_x`, `h_z` — the conditioning (``x``) and transformed-response (``z``) bandwidths (both ``> 0``), supplied.

Reference: B. W. Silverman (1986), [*Density Estimation for Statistics and Data Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall.
"""
function conditional_pushforward_density(summary::Histogram2D, g, x::Real, z::Real;
                                         kernel_x = gaussian_kernel,
                                         kernel_z = gaussian_kernel,
                                         degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                         h_x::Real, h_z::Real)
  h_z > 0 || throw(ArgumentError("response bandwidth h_z must be positive (got h_z = $h_z)"))
  return conditional_expectation(summary, y -> kernel_z((g(y) - z) / h_z) / h_z, x;
                                 h = h_x, kernel = kernel_x, degree = degree)
end

function conditional_pushforward_density(summary::Histogram2D, g, x::Real, zs::AbstractVector;
                                         kernel_x = gaussian_kernel,
                                         kernel_z = gaussian_kernel,
                                         degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                         h_x::Real, h_z::Real)
  h_z > 0 || throw(ArgumentError("response bandwidth h_z must be positive (got h_z = $h_z)"))
  return [conditional_expectation(summary, y -> kernel_z((g(y) - z) / h_z) / h_z, x;
                                  h = h_x, kernel = kernel_x, degree = degree) for z in zs]
end


function _check_cpf_se(degree, h_x, h_z)
  (h_x > 0 && h_z > 0) ||
    throw(ArgumentError("bandwidths h_x, h_z must be positive (got h_x = $h_x, h_z = $h_z)"))
  _check_degree(degree)
end


# The local-fit leverage ℓ(x) at the conditioning point (z-independent, so formed once per SE call).
function _cpf_leverage(summary::Histogram2D, x::Real, kernel_x, degree, h_x)
  c_x, N, _ = _column_response(summary, identity)   # c_x, N do not depend on the response
  return _leverage(c_x, N, x, kernel_x, degree, h_x)
end


# The conditional pushforward density f̂_{g(Y)|X}(z | x) = 𝔼̂[K_{h_z}(g(Y) − z) | X = x] and its standard
# error. f̂_{g(Y)|X} is the local fit of the response K_{h_z}(g(Y) − z), so its variance is the local-fit
# leverage ℓ(x) times the variance of that response, ℓ(x)·(R(K_z) f̂_{g(Y)|X}/h_z − f̂_{g(Y)|X}²) — the
# analytic R(K_z), since the order-2 y-quadrature cannot resolve K_{h_z}². Both terms carry ℓ(x); the
# marginal −f̂²/n has weight 1/n because every sample counts there, which is not so here. At degree 0 this
# reduces asymptotically to the `conditional_density` standard error. Both NaN where the fit is degenerate.
function _cpf_value_and_se(summary::Histogram2D, g, x::Real, z::Real, kernel_x, kernel_z, degree, h_x, h_z, lev)
  f = conditional_pushforward_density(summary, g, x, z; h_x = h_x, h_z = h_z, kernel_x = kernel_x,
                                      kernel_z = kernel_z, degree = degree)
  (isfinite(f) && isfinite(lev)) || return (f, oftype(f, NaN))
  variance = lev * (_quad(u -> kernel_z(u)^2) * f / h_z - f^2)
  return (f, sqrt(max(variance, zero(variance))))
end


"""
    conditional_pushforward_density_standard_error(summary::Histogram2D, g, x, z;
        kernel_x = gaussian_kernel, kernel_z = gaussian_kernel,
        degree = kernel_x isa AsymmetricKernel ? 0 : 1, h_x, h_z)

**Estimand.** The estimand is the pointwise standard error of the conditional pushforward density
`conditional_pushforward_density` ``\\equiv f_{g(Y)\\mid X}``, the sampling error of the density of the transformed response ``g(Y)`` at ``z`` given ``x``.

**Estimator.** The band is the two-term (kernel term minus squared-density term) ``\\sqrt{\\hat\\ell(x) [R(K_z)\\hat f_{g(Y)\\mid X}/h_z - \\hat f_{g(Y)\\mid X}^2]}``:
the local-regression leverage ``\\hat\\ell(x)`` of the fit in ``x`` (correct at any `degree`) times the variance of the
response ``K_{h_z}(g(Y) - z)`` — ``R(K_z) = \\int K_z^2`` by quadrature, the analytic form since the
order-2 ``y``-quadrature cannot resolve the squared narrow kernel ``K_{h_z}^2``. At `degree` 0 it reduces asymptotically to the
[`conditional_density_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_standard_error)
of the transformed response. Only a symmetric `kernel_z` is supported.

Arguments: as for [`conditional_pushforward_density`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_pushforward_density).

"""
function conditional_pushforward_density_standard_error(summary::Histogram2D, g, x::Real, z::Real;
                                                        kernel_x = gaussian_kernel,
                                                        kernel_z = gaussian_kernel,
                                                        degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                                        h_x::Real, h_z::Real)
  _check_cpf_se(degree, h_x, h_z)
  lev = _cpf_leverage(summary, x, kernel_x, degree, h_x)
  return _cpf_value_and_se(summary, g, x, z, kernel_x, kernel_z, degree, h_x, h_z, lev)[2]
end

function conditional_pushforward_density_standard_error(summary::Histogram2D, g, x::Real,
                                                        zs::AbstractVector;
                                                        kernel_x = gaussian_kernel,
                                                        kernel_z = gaussian_kernel,
                                                        degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                                        h_x::Real, h_z::Real)
  _check_cpf_se(degree, h_x, h_z)
  lev = _cpf_leverage(summary, x, kernel_x, degree, h_x)
  return [_cpf_value_and_se(summary, g, x, z, kernel_x, kernel_z, degree, h_x, h_z, lev)[2] for z in zs]
end


"""
    conditional_pushforward_density_with_standard_error(summary::Histogram2D, g, x, z;
        kernel_x = gaussian_kernel, kernel_z = gaussian_kernel,
        degree = kernel_x isa AsymmetricKernel ? 0 : 1, h_x, h_z)

Returns the conditional pushforward density `conditional_pushforward_density` and its standard error from
one fit as the `NamedTuple` ``(\\texttt{value}, \\texttt{standard\\_error})``, sharing the fit's leverage with
[`conditional_pushforward_density_standard_error`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_pushforward_density_standard_error). For a vector
`z` it returns a `NamedTuple` of vectors.

Arguments: as for [`conditional_pushforward_density`](@ref
GBGStreamingStatistics.ConditionalExpectations.conditional_pushforward_density).
"""
function conditional_pushforward_density_with_standard_error(summary::Histogram2D, g, x::Real, z::Real;
                                                             kernel_x = gaussian_kernel,
                                                             kernel_z = gaussian_kernel,
                                                             degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                                             h_x::Real, h_z::Real)
  _check_cpf_se(degree, h_x, h_z)
  lev = _cpf_leverage(summary, x, kernel_x, degree, h_x)
  v, se = _cpf_value_and_se(summary, g, x, z, kernel_x, kernel_z, degree, h_x, h_z, lev)
  return (value = v, standard_error = se)
end

function conditional_pushforward_density_with_standard_error(summary::Histogram2D, g, x::Real,
                                                             zs::AbstractVector;
                                                             kernel_x = gaussian_kernel,
                                                             kernel_z = gaussian_kernel,
                                                             degree::Integer = (kernel_x isa AsymmetricKernel ? 0 : 1),
                                                             h_x::Real, h_z::Real)
  _check_cpf_se(degree, h_x, h_z)
  lev = _cpf_leverage(summary, x, kernel_x, degree, h_x)
  vses = [_cpf_value_and_se(summary, g, x, z, kernel_x, kernel_z, degree, h_x, h_z, lev) for z in zs]
  return (value = [vs[1] for vs in vses], standard_error = [vs[2] for vs in vses])
end

end # module ConditionalExpectations
