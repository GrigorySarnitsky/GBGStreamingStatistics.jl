# ConditionalDensity — §4.3: the conditional density f_{Y|X} and the two-time transition density
# p_τ(y|x), from a joint histogram by the double-kernel estimator (`docs/ARCHITECTURE.md` §4.3).

module ConditionalDensity

using ..Histograms: Histogram2D, _xbincenter, _ybincenter
using ..Kernels: gaussian_kernel, _kernel_weight, AsymmetricKernel

export conditional_density, transition_density


_check_bandwidths(hx::Real, hy::Real) =
  (hx > 0 && hy > 0) ||
    throw(ArgumentError("bandwidths hx, hy must be positive (got hx = $hx, hy = $hy)"))


# Reduce the joint histogram at conditioning value `x` to a 1-D X-kernel-weighted Y-profile:
# `P[k] = Σ_j n_{jk} K_hx(c^x_j, x)` (the gx-weighted count in Y-row k) and `D = Σ_k P[k] = Σ_j N_j
# K_hx(c^x_j, x)` (the marginal-X KDE numerator). The conditional distribution over Y-cells is P/D;
# smoothing it in y gives the density. The per-column weight is the shared `_kernel_weight` atom (symmetric
# or asymmetric `kernelx` alike); the x-weighting cancels in the f_{X,Y}/f_X ratio, so its normalization is
# irrelevant here. Host-side (`Array` is a no-op for CPU summaries).
function _y_profile(summary::Histogram2D, x::Real, hx::Real, kernelx)
  counts = Array(summary.counts)
  nx = summary.nx
  ny = summary.ny
  gx = Vector{Float64}(undef, nx)
  @inbounds for j in 1:nx
    gx[j] = _kernel_weight(_xbincenter(summary, j), x, hx, kernelx)
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
# weighting (it appears in both the profile P and its total D); only the y-direction normalization
# survives, carried by the shared `_kernel_weight` atom (its 1/hy folded in). For a symmetric `kernely`,
# ∫K = 1 makes ∫ f̂(y∣x) dy = 1 exactly; an asymmetric `kernely` (a normalized density in its data argument)
# gives ∫ f̂ dy ≈ 1 — exact only after `_response_normalizer` (as for the marginal `density`).
function _density_at(summary::Histogram2D, P, D, y::Real, hy::Real, kernely)
  D > 0 || return oftype(D, NaN)
  acc = 0.0
  @inbounds for k in 1:summary.ny
    Pk = P[k]
    iszero(Pk) && continue
    acc += Pk * _kernel_weight(_ybincenter(summary, k), y, hy, kernely)
  end
  return acc / D
end


# ∫ f̂(y∣x) dy over the summary's Y-range ∩ the response kernel's support, by the trapezoidal rule on a
# fine grid — the per-x renormalization that makes ∫ f̂ dy = 1 exactly for a non-self-normalizing
# asymmetric response kernel (a symmetric kernel already gives 1, so `normalize` is a no-op there).
function _response_normalizer(summary::Histogram2D, P, D, hy::Real, kernely)
  lower, upper = float(summary.ylower), float(summary.yupper)
  if kernely isa AsymmetricKernel
    lower = max(lower, kernely.support[1])
    upper = min(upper, kernely.support[2])
  end
  npts = 1001
  grid = range(lower, upper; length = npts)
  acc = 0.0
  @inbounds for (i, yy) in enumerate(grid)
    w = (i == 1 || i == npts) ? 0.5 : 1.0
    acc += w * _density_at(summary, P, D, yy, hy, kernely)
  end
  return acc * step(grid)
end


"""
    conditional_density(summary::Histogram2D, x, y;
                        hx, hy, kernelx = gaussian_kernel, kernely = gaussian_kernel,
                        normalize = false)   # ≡ f_{Y|X}

**Estimand.** the conditional density
``f_{Y\\mid X}(y \\mid x) = f_{X,Y}(x, y) / f_X(x)``.

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996). From the joint
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) `summary`, with a conditioning kernel
``K_1 =`` `kernelx` (bandwidth ``h_1 =`` `hx`, in ``x``) and a response kernel ``K_2 =`` `kernely`
(bandwidth ``h_2 =`` `hy`, in ``y``); both ``\\int K = 1``, default `gaussian_kernel`,

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

The ``x``-kernel weighting cancels in the ratio (it appears in both the profile and its total), so the
choice of `kernelx` only weights the local conditioning window. For a **symmetric** response `kernely`
the surviving ``1/h_2`` makes ``\\int \\hat f(y \\mid x)\\,\\mathrm{d}y = 1`` exactly (since ``\\int K =
1``). For a **bounded / positive response** — a droplet radius as ``Y``, say — pass an **asymmetric**
`kernely` ([`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) for ``[0,\\infty)``,
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) for ``[0,1]``), the same boundary-bias
fix the marginal [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) uses: it is a density in
the data argument ``K_2(c^y_k;\\, y, h_2)`` (no mass leaks past the edge) but only ``\\approx`` normalized
in ``y``, so set `normalize = true` to enforce ``\\int \\hat f\\,\\mathrm{d}y = 1`` exactly. (A bounded
*conditioning* variable needs no asymmetric `kernelx`: the ``x``-kernel cancels in the ``f_{X,Y}/f_X``
ratio, so its boundary leakage does not bias ``f_{Y\\mid X}``.) The
conditional MISE rate is slower than a marginal density (``O(n^{-2/3})``). A vector `y` returns the
density curve over those points at the fixed conditioning value `x`; where no data is local to `x` the
estimate is `NaN`.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value (a scalar).
- `y` — the response point(s) at which the density is evaluated (a scalar, or a vector for the curve).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``), supplied.
- `kernelx` — the conditioning-direction (``x``) kernel (any ``\\int K = 1`` callable; default
  `gaussian_kernel`). Its normalization cancels in the ``f_{X,Y}/f_X`` ratio, so a symmetric kernel is
  almost always the right choice.
- `kernely` — the response-direction (``y``) kernel (default `gaussian_kernel`). A **symmetric** kernel
  gives ``\\int\\hat f\\,\\mathrm dy = 1`` exactly; an **asymmetric**
  [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel)
  (`gamma_kernel` / `beta_kernel`), for a bounded / positive response, is boundary-bias-free but needs
  `normalize = true` for an exact integral.
- `normalize` — when `true`, divide by ``\\int\\hat f(y\\mid x)\\,\\mathrm dy`` over the response range so
  the curve integrates to 1 **exactly**: the finishing step for an asymmetric `kernely`, a no-op for a
  symmetric one.

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
function conditional_density(summary::Histogram2D, x::Real, y::Real; hx::Real, hy::Real,
                             kernelx = gaussian_kernel, kernely = gaussian_kernel,
                             normalize::Bool = false)
  _check_bandwidths(hx, hy)
  P, D = _y_profile(summary, x, hx, kernelx)
  f = _density_at(summary, P, D, y, hy, kernely)
  (normalize && D > 0) || return f
  return f / _response_normalizer(summary, P, D, hy, kernely)
end

function conditional_density(summary::Histogram2D, x::Real, ys::AbstractVector; hx::Real, hy::Real,
                             kernelx = gaussian_kernel, kernely = gaussian_kernel,
                             normalize::Bool = false)
  _check_bandwidths(hx, hy)
  P, D = _y_profile(summary, x, hx, kernelx)
  Z = (normalize && D > 0) ? _response_normalizer(summary, P, D, hy, kernely) : 1.0
  return [_density_at(summary, P, D, y, hy, kernely) / Z for y in ys]
end


"""
    transition_density(summary::Histogram2D, x, y;
                       hx, hy, kernelx = gaussian_kernel, kernely = gaussian_kernel,
                       normalize = false)   # ≡ p_τ(y|x)

**Estimand.** the **two-time transition density**
``p_\\tau(y \\mid x) = f_{X_{t+\\tau} \\mid X_t}(y \\mid x)`` — the conditional density of the state a
lag ``\\tau`` later, given the present state (nonstationary: ``p_{t,\\tau}``). For the
Ornstein–Uhlenbeck process this is Gaussian in closed form.

**Estimator.** the [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)
double-kernel estimator applied to the per-lag joint histogram of pairs
``\\bigl(X_t,\\, X_{t+\\tau}\\bigr)`` (binned by
[`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!) at the summary's lag
``\\tau``); the computation is identical — the name marks the two-time reading. The bandwidths `hx`
(conditioning, in ``x``), `hy` (response, in ``y``), the kernels `kernelx` / `kernely`, and `normalize`
are as for
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density). A vector
`y` returns the transition-density curve at the fixed present state `x`.
"""
transition_density(summary::Histogram2D, x::Real, y; hx::Real, hy::Real,
                   kernelx = gaussian_kernel, kernely = gaussian_kernel, normalize::Bool = false) =
  conditional_density(summary, x, y; hx = hx, hy = hy, kernelx = kernelx, kernely = kernely,
                      normalize = normalize)

end # module ConditionalDensity
