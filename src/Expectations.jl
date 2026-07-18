# Expectations — §4.B: the power-sum quadrature readouts, and the recovered marginal moments they
# generalize. A smooth functional of X read straight off the stored power-sums, with no kernel and no
# query point: the recovered moments `mean` / `var` / `std` / `central_moment` (the g = x^p case), the
# expectation E[g(X)] of a general smooth g, its sampling standard error, and the pushforward density of
# Y = g(X) (`docs/ARCHITECTURE.md` §4.B). All route through the one shared contraction atom `_quadrature`
# — the moments supplying closed-form (binomial) coefficients, the general `expectation` a local
# degree-`order` polynomial interpolant of g (dependency-free; exact on the monomials, so it reduces to
# the shift formula). Pure host-side numerics over the recovered summary.

module Expectations

import Statistics: mean, var, std

using LinearAlgebra: inv, mul!

using ..Histograms: Histogram1D, _bincenter, _order, _total
using ..Kernels: gaussian_kernel


export mean, var, std, central_moment,
       expectation, expectation_standard_error, expectation_with_standard_error,
       pushforward_density, pushforward_density_standard_error, pushforward_density_with_standard_error


# --- Power-sum quadrature: the shared atom for every global (non-local) readout ---
# E[g(X)] = (1/n) Σ_j Σ_{p=0}^{order} a_p(c_j) S_p^(j), with S_0^(j) = n_j and a_p(c_j) the local
# coefficients of the test function's degree-≤order polynomial model about the bin center c_j (a[1] the
# constant, weighting the count; a[p+1] multiplying (x−c_j)^p). `coeffs!(a, c_j)` fills the reused
# length-(order+1) buffer; it is called only on non-empty bins, so g is never probed off the data. The
# recovered moments are the g = x^p / (x−μ)^p case (closed-form binomial coefficients — the shift
# formula); the general `expectation` supplies interpolation coefficients for a general smooth g. Reads
# the summary to the host first (a no-op for CPU arrays). Returns E[g] (NaN for an empty summary).
function _quadrature(h::Histogram1D, coeffs!)
  P = _order(h)
  counts = Array(h.counts)
  S = Array(h.powersums)
  a = zeros(Float64, P + 1)          # reused per-bin coefficient buffer (a[p+1] ↔ (x−c_j)^p)
  n = 0.0
  acc = 0.0
  @inbounds for j in eachindex(counts)
    nj = counts[j]
    iszero(nj) && continue
    n += nj
    coeffs!(a, _bincenter(h, j))
    acc += a[1] * nj                 # p = 0 term: S_0^(j) = n_j
    for p in 1:P
      acc += a[p + 1] * S[p, j]
    end
  end
  return iszero(n) ? NaN : acc / n
end


# Local coefficients of the raw monomial x^m about c_j: (c_j + u)^m = Σ_k C(m,k) c_j^{m-k} u^k, so
# a[k+1] = C(m,k) c_j^{m-k} for k ≤ m (0 above). The g = x^m instance of the quadrature.
function _monomial_coeffs!(a, cj, m, P)
  fill!(a, 0.0)
  @inbounds for k in 0:min(m, P)
    a[k + 1] = binomial(m, k) * cj^(m - k)
  end
  return a
end


# Local coefficients of the central monomial (x−μ)^p about c_j: (c_j − μ + u)^p, so
# a[k+1] = C(p,k)(c_j − μ)^{p-k}. The g = (x−μ)^p instance (the shift formula for central moments).
function _central_coeffs!(a, cj, mu, p, P)
  fill!(a, 0.0)
  d = cj - mu
  @inbounds for k in 0:min(p, P)
    a[k + 1] = binomial(p, k) * d^(p - k)
  end
  return a
end


_sample_mean(h::Histogram1D) =
  (P = _order(h); _quadrature(h, (a, cj) -> _monomial_coeffs!(a, cj, 1, P)))

function _sample_central_moment(h::Histogram1D, p::Integer)
  2 ≤ p ≤ _order(h) ||
    throw(ArgumentError("central-moment order p must satisfy 2 ≤ p ≤ $(_order(h)) (got $p)"))
  mu = _sample_mean(h)
  isnan(mu) && return oftype(mu, NaN)
  P = _order(h)
  return _quadrature(h, (a, cj) -> _central_coeffs!(a, cj, mu, p, P))
end

_sample_variance(h::Histogram1D) = _sample_central_moment(h, 2)
_sample_std(h::Histogram1D) = sqrt(_sample_variance(h))


"""
    mean(summary::Histogram1D)

**Estimand.** The population mean ``\\mathbb{E}[X]``.

**Estimator.** The sample mean ``\\bar X = \\frac{1}{n} \\sum_i X_i`` — recovered **exactly** (not a smoothed
estimate) host-side from the per-bin counts and first power-sums by the shift formula (the power-sums use
the actual ``X_i``, not the bin center; the ``g = x`` case of [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation)). Extends `Statistics.mean`; an empty summary gives `NaN`.
"""
mean(summary::Histogram1D) = _sample_mean(summary)


"""
    var(summary::Histogram1D; corrected = true)

**Estimand.** The population variance ``\\operatorname{Var}(X) = \\mathbb{E}\\bigl[(X - \\mathbb{E}[X])^2\\bigr]``.

**Estimator.** The recovered sample variance (``\\bar X`` the recovered mean). With `corrected = true`
(the default, matching `Statistics.var`) the **unbiased** ``\\frac{1}{n-1} \\sum_i (X_i - \\bar X)^2``
(Bessel's correction; `NaN` for ``n \\le 1``); `corrected = false` the ``\\frac{1}{n} \\sum_i (X_i - \\bar X)^2``
order-2 central moment (`central_moment` at ``p = 2``). Exact via the shift formula. Extends
`Statistics.var`.
"""
function var(summary::Histogram1D; corrected::Bool = true)
  m2 = _sample_variance(summary)
  corrected || return m2
  n = _total(summary)
  return n > 1 ? m2 * n / (n - 1) : oftype(m2, NaN)
end


"""
    std(summary::Histogram1D; corrected = true)

**Estimand.** The population standard deviation ``\\sigma = \\sqrt{\\operatorname{Var}(X)}``.

**Estimator.** The recovered sample standard deviation ``\\sqrt{s^2}``, where ``s^2`` is the recovered
variance (`var`); the `corrected` keyword has the same meaning (the default `true` divides by ``n - 1``).
Exact via the shift formula. Extends `Statistics.std`.
"""
std(summary::Histogram1D; corrected::Bool = true) = sqrt(var(summary; corrected = corrected))


"""
    central_moment(summary::Histogram1D, p)

**Estimand.** The population ``p``-th central moment ``\\mu_p = \\mathbb{E}\\bigl[(X - \\mathbb{E}[X])^p\\bigr]``.

**Estimator.** The sample central moment ``\\frac{1}{n} \\sum_i (X_i - \\bar X)^p`` (``\\bar X`` the recovered
mean), recovered **exactly** from the power-sums by the shift formula — not a smoothed estimate — for
``2 \\le p \\le \\texttt{order}``. ``p = 2`` is the uncorrected variance (`var` with `corrected = false`);
even orders recover stably, odd orders (``\\approx 0``) to weaker relative precision. This is the
``g = (x - \\bar X)^p`` case of the power-sum quadrature [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation), with which it shares the contraction atom.
"""
central_moment(summary::Histogram1D, p::Integer) = _sample_central_moment(summary, p)


# Chebyshev–Lobatto nodes on [-1, 1] and the inverse of their Vandermonde matrix, for the degree-`P`
# local-polynomial model of g on a bin. Uniform bins ⇒ one matrix serves every bin (the nodes are the same
# fixed offsets, scaled by the half-width). Interpolating in the scaled variable s = (x − c_j)/(Δ/2) keeps
# the Vandermonde Δ-independent and well-conditioned at the low orders used here. Called once per readout.
function _interp_matrix(P::Integer)
  nodes = [cos(π * (m - 1) / P) for m in 1:P + 1]              # P+1 Chebyshev–Lobatto points on [-1, 1]
  vandermonde = [nodes[m]^(k - 1) for m in 1:P + 1, k in 1:P + 1]
  return nodes, inv(vandermonde)
end


"""
    expectation(summary::Histogram1D, g)   # ≡ 𝔼[g(X)]

**Estimand.** The population expectation of a smooth functional of ``X``,
``\\mathbb{E}[g(X)] = \\int g(x)  f_X(x)\\,\\mathrm{d}x``.

**Estimator** (power-sum quadrature). From the [`Histogram1D`](@ref
GBGStreamingStatistics.Histograms.Histogram1D) `summary`, whose bin ``j`` (center ``c_j``) carries the
count ``n_j`` and its per-bin power-sums ``S_p^{(j)}`` about the center (with ``S_0^{(j)} = n_j``),
Taylor-expanding ``g`` about each center gives

```math
\\widehat{\\mathbb{E}}[g(X)] = \\frac{1}{n} \\sum_j \\sum_{p=0}^{\\texttt{order}}
  \\frac{g^{(p)}(c_j)}{p!}  S_p^{(j)} ,
```

which is **exact** for a polynomial ``g`` of degree ``\\le \\texttt{order}`` (the remainder terminates —
this is the shift formula behind `mean` / `var` / `central_moment`, the ``g = x^p`` case) and
``O(\\Delta^{ \\texttt{order}+1})`` in the bin width ``\\Delta`` for a general smooth ``g``. The local
coefficients ``g^{(p)}(c_j)/p!`` are read off a degree-``\\texttt{order}`` polynomial interpolant of ``g``
on each bin (no automatic differentiation), so any callable `g` works; `g` is probed only on non-empty
bins, hence only within the data range. The estimand needs ``g`` smooth on the bin scale and the mass
inside the summary's range. An empty summary gives `NaN`.

Arguments:

- `summary` — the [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) of ``X``.
- `g` — the functional to average: any callable ``x \\mapsto g(x)`` (smooth on the bin scale). The moments
  are the ``g = x^p`` case (`mean` is `g = identity`); `g = sin`, `g = exp`, etc. are the general case.

```jldoctest
julia> h = Histogram1D(-4.0, 4.0, 800);

julia> observe!(h, [-1.0, 0.0, 2.0]);   # the g = identity case is the recovered mean E[X] = 1/3

julia> expectation(h, identity) ≈ mean(h)
true
```
"""
function expectation(summary::Histogram1D, g)
  P = _order(summary)
  halfwidth = 0.5 / summary.scale                             # Δ/2 = 0.5 (upper − lower) / nbins
  nodes, vandermonde_inv = _interp_matrix(P)
  gsamples = Vector{Float64}(undef, P + 1)                    # reused g-sample buffer
  scaled = Vector{Float64}(undef, P + 1)                      # reused coefficient buffer (in s)
  invpow = [halfwidth^(-p) for p in 0:P]                      # halfwidth^{-p}: converts s-coeffs → u-coeffs
  coeffs! = function (a, cj)
    @inbounds for m in 1:P + 1
      gsamples[m] = g(cj + nodes[m] * halfwidth)
    end
    mul!(scaled, vandermonde_inv, gsamples)                  # coefficients in s = (x − c_j)/halfwidth
    @inbounds for p in 0:P
      a[p + 1] = scaled[p + 1] * invpow[p + 1]               # a_p = b_p / halfwidth^p (s → u = x − c_j)
    end
    return a
  end
  return _quadrature(summary, coeffs!)
end


"""
    expectation_standard_error(summary::Histogram1D, g)

**Estimand.** The sampling standard error of the estimator [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation),
``\\mathrm{SE}\\bigl(\\widehat{\\mathbb{E}}[g(X)]\\bigr) = \\sqrt{\\operatorname{Var}(g(X))/n}`` with
``\\operatorname{Var}(g(X)) = \\mathbb{E}[g(X)^2] - \\mathbb{E}[g(X)]^2``.

**Estimator.** ``\\sqrt{(\\widehat{\\mathbb{E}}[g^2] - \\widehat{\\mathbb{E}}[g]^2)/n}`` — both moments the same
power-sum quadrature (``\\mathbb{E}[g^2]`` is [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation) with the smooth integrand ``g^2``), ``n`` the sample
count. This is the plug-in sampling error of the empirical mean of ``g(X)``; it does not include the
``O(\\Delta^{ \\texttt{order}+1})`` quadrature bias, which is separate. `NaN` for an empty summary.

Arguments: as for [`expectation`](@ref GBGStreamingStatistics.Expectations.expectation) — `summary`, `g`.
"""
function expectation_standard_error(summary::Histogram1D, g)
  n = _total(summary)
  n > 0 || return NaN
  eg = expectation(summary, g)
  eg2 = expectation(summary, x -> g(x)^2)
  variance = eg2 - eg^2
  return sqrt(max(variance, zero(variance)) / n)
end


"""
    expectation_with_standard_error(summary::Histogram1D, g)

The estimate and its standard error from one power-sum quadrature pass — the `NamedTuple`
``(\\texttt{value}, \\texttt{standard\\_error})``, with `value` the [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation) ``\\widehat{\\mathbb{E}}[g(X)]`` and `standard_error` its
[`expectation_standard_error`](@ref GBGStreamingStatistics.Expectations.expectation_standard_error). They
describe the same estimator by construction, and ``\\widehat{\\mathbb{E}}[g]`` is computed once and shared.

Arguments: as for [`expectation`](@ref GBGStreamingStatistics.Expectations.expectation).

```jldoctest
julia> using GBGStreamingStatistics

julia> h = Histogram1D(-8.0, 8.0, 800); observe!(h, [-1.0, 0.0, 0.5, 1.0, 2.0]);

julia> r = expectation_with_standard_error(h, cos); (r.value isa Float64, r.standard_error > 0)
(true, true)
```
"""
function expectation_with_standard_error(summary::Histogram1D, g)
  n = _total(summary)
  eg = expectation(summary, g)
  n > 0 || return (value = eg, standard_error = NaN)
  eg2 = expectation(summary, x -> g(x)^2)
  variance = eg2 - eg^2
  return (value = eg, standard_error = sqrt(max(variance, zero(variance)) / n))
end


"""
    pushforward_density(summary::Histogram1D, g, y;
                        h, kernel = gaussian_kernel)   # ≡ f_{g(X)}

**Estimand.** The density of the transformed variable ``Y = g(X)``, the **pushforward**
``f_{g(X)}(y) = \\sum_{x: g(x)=y} f_X(x) / |g'(x)|`` (the change-of-variables sum over the preimages
of ``y``).

**Estimator.** The kernel density estimate of the pushforward, written as an expectation and evaluated by
[`expectation`](@ref GBGStreamingStatistics.Expectations.expectation):

```math
\\hat f_{g(X)}(y) = \\frac{1}{n} \\sum_i K_h(g(X_i) - y)
  = \\widehat{\\mathbb{E}}\\bigl[  K_h(g(X) - y)  \\bigr] ,
```

so the composed test function ``x \\mapsto K_h(g(x) - y)`` is smoothed by the power-sum quadrature — no
inverse ``g^{-1}`` and no explicit Jacobian are formed; a **non-monotone** ``g`` is handled automatically
(every preimage contributes), and a flat spot ``g'(x) = 0`` shows up as the expected (integrable) caustic
pile-up. For a symmetric `kernel` (``\\int K = 1``) the surviving ``1/h`` makes ``\\int \\hat f\\,\\mathrm
dy = 1`` exactly. Needs ``g`` smooth on the bin scale and the response bandwidth ``h`` at least of order
``|g'| \\Delta`` (so ``K_h \\circ g`` is resolved across a bin). A vector `y` returns the density curve
over those points.

Arguments:

- `summary` — the [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) of ``X``.
- `g` — the transform defining ``Y = g(X)`` (any callable, smooth on the bin scale).
- `y` — the response point(s) at which the pushforward density is evaluated (a scalar, or a vector for the
  curve).
- `h` — the response-direction bandwidth (``h > 0``), supplied.
- `kernel` — the response-direction kernel (any ``\\int K = 1`` callable; default
  [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)).

Reference: B. W. Silverman (1986), [*Density Estimation for Statistics and Data Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall.
"""
function pushforward_density(summary::Histogram1D, g, y::Real; h::Real, kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  return expectation(summary, x -> kernel((g(x) - y) / h) / h)
end

function pushforward_density(summary::Histogram1D, g, ys::AbstractVector;
                             h::Real, kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  return [expectation(summary, x -> kernel((g(x) - y) / h) / h) for y in ys]
end


"""
    pushforward_density_standard_error(summary::Histogram1D, g, y;
                                       h, kernel = gaussian_kernel)

**Estimand.** The pointwise standard error of the pushforward density [`pushforward_density`](@ref
GBGStreamingStatistics.Expectations.pushforward_density) ``\\equiv f_{g(X)}``,
``\\mathrm{SE}(y) = \\sqrt{\\operatorname{Var}\\hat f_{g(X)}(y)}`` with the two-term kernel-density variance
``\\operatorname{Var}\\hat f_{g(X)}(y) \\approx R(K) f_{g(X)}(y)/(nh) - f_{g(X)}(y)^2/n`` (``R(K) = \\int K^2``).

**Estimator.** Because ``\\hat f_{g(X)}(y) = \\widehat{\\mathbb{E}}[K_h(g(X) - y)]`` is a sample mean of the
composed kernel, its standard error is exactly the
[`expectation_standard_error`](@ref GBGStreamingStatistics.Expectations.expectation_standard_error) of that
integrand — ``\\sqrt{(\\widehat{\\mathbb{E}}[K_h(g(X)-y)^2] - \\hat f_{g(X)}(y)^2)/n}``, which reduces to the
two-term form above.

Arguments: as for [`pushforward_density`](@ref GBGStreamingStatistics.Expectations.pushforward_density) —
`summary`, `g`, `y`; `h`, `kernel`.
"""
function pushforward_density_standard_error(summary::Histogram1D, g, y::Real; h::Real,
                                            kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  return expectation_standard_error(summary, x -> kernel((g(x) - y) / h) / h)
end

function pushforward_density_standard_error(summary::Histogram1D, g, ys::AbstractVector; h::Real,
                                            kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  return [expectation_standard_error(summary, x -> kernel((g(x) - y) / h) / h) for y in ys]
end


"""
    pushforward_density_with_standard_error(summary::Histogram1D, g, y;
                                            h, kernel = gaussian_kernel)

[`pushforward_density`](@ref GBGStreamingStatistics.Expectations.pushforward_density) and its standard
error from one quadrature pass — the `NamedTuple` ``(\\texttt{value}, \\texttt{standard\\_error})``, sharing
the composed-kernel integrand with [`pushforward_density_standard_error`](@ref
GBGStreamingStatistics.Expectations.pushforward_density_standard_error). For a vector `y`, a `NamedTuple` of
vectors.

Arguments: as for [`pushforward_density`](@ref GBGStreamingStatistics.Expectations.pushforward_density).
"""
function pushforward_density_with_standard_error(summary::Histogram1D, g, y::Real; h::Real,
                                                 kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  return expectation_with_standard_error(summary, x -> kernel((g(x) - y) / h) / h)
end

function pushforward_density_with_standard_error(summary::Histogram1D, g, ys::AbstractVector; h::Real,
                                                 kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  rs = [expectation_with_standard_error(summary, x -> kernel((g(x) - y) / h) / h) for y in ys]
  return (value = [r.value for r in rs], standard_error = [r.standard_error for r in rs])
end


end # module Expectations
