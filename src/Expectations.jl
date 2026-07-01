# Expectations — §4.B: the power-sum quadrature readouts. A smooth functional of X read straight off the
# stored power-sums, with no kernel and no query point: the expectation E[g(X)] of a general smooth g, its
# sampling standard error, and the pushforward density of Y = g(X) (`docs/ARCHITECTURE.md` §4.B). All three
# route through the shared `Histograms._quadrature` atom — the same contraction the recovered moments use
# (the moments are the g = x^p case). A general g's per-bin coefficients come from a local degree-`order`
# polynomial interpolant (dependency-free; exact for polynomials of degree ≤ order, so it reduces to the
# shift formula on the monomials). Pure host-side numerics over the recovered summary.

module Expectations

using LinearAlgebra: inv, mul!

using ..Histograms: Histogram1D, _quadrature, _order, _total
using ..Kernels: gaussian_kernel


export expectation, expectation_standard_error, pushforward_density


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

**Estimand.** the population expectation of a smooth functional of ``X``,
``\\mathbb{E}[g(X)] = \\int g(x)\\, f_X(x)\\,\\mathrm{d}x``.

**Estimator** (power-sum quadrature). From the [`Histogram1D`](@ref
GBGStreamingStatistics.Histograms.Histogram1D) `summary`, whose bin ``j`` (center ``c_j``) carries the
count ``n_j`` and the power-sums ``S_p^{(j)} = \\sum_{i\\in\\text{bin}}(X_i - c_j)^p`` about the center
(``S_0^{(j)} = n_j``), Taylor-expanding ``g`` about each center gives

```math
\\hat{\\mathbb{E}}[g(X)] \\;=\\; \\frac1n \\sum_j \\sum_{p=0}^{\\texttt{order}}
  \\frac{g^{(p)}(c_j)}{p!}\\, S_p^{(j)} ,
```

which is **exact** for a polynomial ``g`` of degree ``\\le \\texttt{order}`` (the remainder terminates —
this is the shift formula behind `mean` / `var` / `central_moment`, the ``g = x^p`` case) and
``O(\\Delta^{\\,\\texttt{order}+1})`` in the bin width ``\\Delta`` for a general smooth ``g``. The local
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

**Estimand.** the sampling standard error of the estimator [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation),
``\\mathrm{SE}\\bigl(\\hat{\\mathbb{E}}[g(X)]\\bigr) = \\sqrt{\\operatorname{Var}(g(X))/n}`` with
``\\operatorname{Var}(g(X)) = \\mathbb{E}[g(X)^2] - \\mathbb{E}[g(X)]^2``.

**Estimator.** ``\\sqrt{(\\hat{\\mathbb{E}}[g^2] - \\hat{\\mathbb{E}}[g]^2)/n}`` — both moments the same
power-sum quadrature (``\\mathbb{E}[g^2]`` is [`expectation`](@ref
GBGStreamingStatistics.Expectations.expectation) with the smooth integrand ``g^2``), ``n`` the sample
count. This is the plug-in sampling error of the empirical mean of ``g(X)``; it does not include the
``O(\\Delta^{\\,\\texttt{order}+1})`` quadrature bias, which is separate. `NaN` for an empty summary.

Arguments:

- `summary` — the [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) of ``X``.
- `g` — the functional whose average's standard error is wanted (as in [`expectation`](@ref
  GBGStreamingStatistics.Expectations.expectation)).
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
    pushforward_density(summary::Histogram1D, g, y;
                        h, kernel = gaussian_kernel)   # ≡ f_{g(X)}

**Estimand.** the density of the transformed variable ``Y = g(X)``, the **pushforward**
``f_{g(X)}(y) = \\sum_{x:\\,g(x)=y} f_X(x)\\,/\\,|g'(x)|`` (the change-of-variables sum over the preimages
of ``y``).

**Estimator.** the kernel density estimate of the pushforward, written as an expectation and evaluated by
[`expectation`](@ref GBGStreamingStatistics.Expectations.expectation):

```math
\\hat f_{g(X)}(y) \\;=\\; \\frac1n \\sum_i K_h\\bigl(g(X_i) - y\\bigr)
  \\;=\\; \\hat{\\mathbb{E}}\\Bigl[\\, K_h\\bigl(g(X) - y\\bigr) \\Bigr] ,
```

so the composed test function ``x \\mapsto K_h(g(x) - y)`` is smoothed by the power-sum quadrature — no
inverse ``g^{-1}`` and no explicit Jacobian are formed; a **non-monotone** ``g`` is handled automatically
(every preimage contributes), and a flat spot ``g'(x) = 0`` shows up as the expected (integrable) caustic
pile-up. For a symmetric `kernel` (``\\int K = 1``) the surviving ``1/h`` makes ``\\int \\hat f\\,\\mathrm
dy = 1`` exactly. Needs ``g`` smooth on the bin scale and the response bandwidth ``h`` at least of order
``|g'|\\,\\Delta`` (so ``K_h \\circ g`` is resolved across a bin). A vector `y` returns the density curve
over those points.

Arguments:

- `summary` — the [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D) of ``X``.
- `g` — the transform defining ``Y = g(X)`` (any callable, smooth on the bin scale).
- `y` — the response point(s) at which the pushforward density is evaluated (a scalar, or a vector for the
  curve).
- `h` — the response-direction bandwidth (``h > 0``), supplied.
- `kernel` — the response-direction kernel (any ``\\int K = 1`` callable; default `gaussian_kernel`).

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


end # module Expectations
