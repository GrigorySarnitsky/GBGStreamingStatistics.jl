# ConditionalDensity — §4.3: the conditional density f_{Y|X} and the two-time transition density
# p_τ(y|x), from a joint histogram by the double-kernel estimator (`docs/ARCHITECTURE.md` §4.3).

module ConditionalDensity

using ..Histograms: Histogram2D, _xbincenter, _ybincenter

export conditional_density, transition_density


@inline _gaussian(u) = exp(-u^2 / 2)   # unnormalised; the √(2π) is restored once, in the y-direction


_check_bandwidths(hx::Real, hy::Real) =
  (hx > 0 && hy > 0) ||
    throw(ArgumentError("bandwidths hx, hy must be positive (got hx = $hx, hy = $hy)"))


# Reduce the joint histogram at conditioning value `x` to a 1-D X-kernel-weighted Y-profile:
# `P[k] = Σ_j n_{jk} K_hx(c^x_j − x)` (the gx-weighted count in Y-row k) and `D = Σ_k P[k] = Σ_j N_j
# K_hx(c^x_j − x)` (the marginal-X KDE numerator). The conditional distribution over Y-cells is then
# `P/D`; smoothing it in y gives the density. Host-side (`Array` is a no-op for CPU summaries).
function _y_profile(summary::Histogram2D, x::Real, hx::Real)
  counts = Array(summary.counts)
  nx = summary.nx
  ny = summary.ny
  invhx = inv(float(hx))
  gx = Vector{Float64}(undef, nx)
  @inbounds for j in 1:nx
    gx[j] = _gaussian((_xbincenter(summary, j) - x) * invhx)
  end
  P = zeros(Float64, ny)
  @inbounds for k in 1:ny
    base = (k - 1) * nx
    acc = 0.0
    for j in 1:nx
      c = counts[base + j]
      iszero(c) && continue
      acc += c * gx[j]
    end
    P[k] = acc
  end
  return P, sum(P)
end


# Smooth the Y-profile to the density at `y`. The double-kernel ratio f_{X,Y}/f_X cancels the X-kernel
# normalisation (the 1/hx and one 1/√(2π)); only the y-direction normalisation 1/(√(2π)·hy) survives,
# which is exactly what makes ∫ f̂(y∣x) dy = 1.
function _density_at(summary::Histogram2D, P, D, y::Real, hy::Real)
  D > 0 || return oftype(D, NaN)
  invhy = inv(float(hy))
  acc = 0.0
  @inbounds for k in 1:summary.ny
    Pk = P[k]
    iszero(Pk) && continue
    acc += Pk * _gaussian((_ybincenter(summary, k) - y) * invhy)
  end
  return acc / (sqrt(2π) * hy * D)
end


"""
    conditional_density(summary::Histogram2D, x, y; hx, hy)   # ≡ f_{Y|X}

**Estimand.** the conditional density
``f_{Y\\mid X}(y \\mid x) = f_{X,Y}(x, y) / f_X(x)``.

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996). From the joint
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) `summary`, with the fixed Gaussian
kernel ``K`` (``\\int K = 1``; there is no `kernel` argument), conditioning bandwidth ``h_1 =`` `hx`
(in ``x``) and response bandwidth ``h_2 =`` `hy` (in ``y``),

```math
\\hat f(y \\mid x) =
  \\frac{\\sum_i K_{h_1}(X_i - x)\\, K_{h_2}(Y_i - y)}{\\sum_i K_{h_1}(X_i - x)} ,
```

evaluated from the bin counts ``n_{jk}`` at cell centers ``(c^x_j, c^y_k)`` (the binned form, with
``N_j = \\sum_k n_{jk}`` the per-column count):

```math
\\hat f(y \\mid x) =
  \\frac{\\sum_{j,k} n_{jk}\\, K_{h_1}(c^x_j - x)\\, K_{h_2}(c^y_k - y)}
       {\\sum_{j} N_j\\, K_{h_1}(c^x_j - x)} .
```

The ``x``-kernel normalization cancels in the ratio, leaving exactly the ``1/(\\sqrt{2\\pi}\\,h_2)``
that makes ``\\int \\hat f(y \\mid x)\\,\\mathrm{d}y = 1``. The conditional MISE rate is slower than a
marginal density (``O(n^{-2/3})``). The bandwidths are supplied (the response-aware local plug-in is
forthcoming, §4.5); a vector `y` returns the density curve over those points at the fixed conditioning
value `x`. Where no data is local to `x` the estimate is `NaN`.

```jldoctest
julia> g = Histogram2D(-1.0, 1.0, 2, -1.0, 1.0, 2);

julia> observe_pairs!(g, [0.4], [0.4]);   # one pair → cell with centers (0.5, 0.5)

julia> round(conditional_density(g, 0.5, 0.5; hx = 0.5, hy = 0.2); digits = 4)
1.9947
```

References:

- M. Rosenblatt (1969), *Conditional probability density and regression estimators*, in [*Multivariate Analysis II*](https://openlibrary.org/works/OL7299345W).
- R. J. Hyndman, D. M. Bashtannyk & G. K. Grunwald (1996), [*Estimating and visualizing conditional densities*](https://doi.org/10.1080/10618600.1996.10474715), J. Comput. Graph. Statist. 5(4).
"""
function conditional_density(summary::Histogram2D, x::Real, y::Real; hx::Real, hy::Real)
  _check_bandwidths(hx, hy)
  P, D = _y_profile(summary, x, hx)
  return _density_at(summary, P, D, y, hy)
end

function conditional_density(summary::Histogram2D, x::Real, ys::AbstractVector; hx::Real, hy::Real)
  _check_bandwidths(hx, hy)
  P, D = _y_profile(summary, x, hx)
  return [_density_at(summary, P, D, y, hy) for y in ys]
end


"""
    transition_density(summary::Histogram2D, x, y; hx, hy)   # ≡ p_τ(y|x)

**Estimand.** the **two-time transition density**
``p_\\tau(y \\mid x) = f_{X_{t+\\tau} \\mid X_t}(y \\mid x)`` — the conditional density of the state a
lag ``\\tau`` later, given the present state (nonstationary: ``p_{t,\\tau}``). For the
Ornstein–Uhlenbeck oracle this is Gaussian in closed form (the natural end-to-end test).

**Estimator.** the [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)
double-kernel estimator applied to the per-lag joint histogram of pairs
``\\bigl(X_t,\\, X_{t+\\tau}\\bigr)`` (binned by
[`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!) at the summary's lag
``\\tau``); the computation is identical — the name marks the two-time reading. A vector `y` returns
the transition-density curve at the fixed present state `x`.
"""
transition_density(summary::Histogram2D, x::Real, y; hx::Real, hy::Real) =
  conditional_density(summary, x, y; hx = hx, hy = hy)

end # module ConditionalDensity
