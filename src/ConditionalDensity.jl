# ConditionalDensity — §4.3: the conditional density f_{Y|X} and the two-time transition density
# p_τ(y|x), from a joint histogram by the double-kernel estimator (`docs/ARCHITECTURE.md` §4.3).

module ConditionalDensity

using ..Histograms: Histogram2D, _xbincenter, _ybincenter, _total
using ..Kernels: gaussian_kernel, _kernel_weight, AsymmetricKernel, _quad

export conditional_density, conditional_density_standard_error, conditional_density_with_standard_error,
       transition_density, transition_density_standard_error, transition_density_with_standard_error


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

**Estimand.** The conditional density
``f_{Y\\mid X}(y \\mid x) = f_{X,Y}(x, y) / f_X(x)``.

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996). From the joint
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) `summary`, with a conditioning kernel
``K_1`` (`kernelx`, of bandwidth ``h_1`` given by `hx`, in ``x``) and a response kernel ``K_2`` (`kernely`,
of bandwidth ``h_2`` given by `hy`, in ``y``); both ``\\int K = 1``, default
[`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel),

```math
\\hat f(y \\mid x) =
  \\frac{\\sum_i K_{h_1}(X_i - x)  K_{h_2}(Y_i - y)}{\\sum_i K_{h_1}(X_i - x)} ,
```

evaluated from the bin counts ``n_{jk}`` at cell centers ``(c^x_j, c^y_k)`` (the binned form, with
``N_j = \\sum_k n_{jk}`` the per-column count):

```math
\\hat f(y \\mid x) =
  \\frac{\\sum_{j,k} n_{jk}  K_{h_1}(c^x_j - x)  K_{h_2}(c^y_k - y)}
       {\\sum_{j} N_j  K_{h_1}(c^x_j - x)} .
```

The ``x``-kernel weighting cancels in the ``f_{X,Y}/f_X`` ratio (it appears in both the profile and its
total), so the choice of `kernelx` only weights the local conditioning window — a bounded *conditioning*
variable thus needs no asymmetric `kernelx`. For a **symmetric** response `kernely`
the surviving ``1/h_2`` makes ``\\int \\hat f(y \\mid x)\\,\\mathrm{d}y = 1`` exactly (since ``\\int K =
1``). For a **bounded / positive response** — a droplet radius as ``Y``, say — pass an **asymmetric**
`kernely` ([`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) for ``[0,\\infty)``,
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) for ``[0,1]``), the same boundary-bias
fix the marginal [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) uses: it is a density in
the data argument ``K_2(c^y_k;  y, h_2)`` (no mass leaks past the edge) but only ``\\approx`` normalized
in ``y``, so set `normalize = true` to enforce ``\\int \\hat f\\,\\mathrm{d}y = 1`` exactly. The
conditional MISE rate is slower than a marginal density (``O(n^{-2/3})``). A vector `y` returns the
density curve over those points at the fixed conditioning value `x`; where no data is local to `x` the
estimate is `NaN`.

# Arguments

- `summary` — the joint `Histogram2D` of ``(X, Y)``.
- `x` — the conditioning value (a scalar).
- `y` — the response point(s) at which the density is evaluated (a scalar, or a vector for the curve).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``), supplied.
- `kernelx` — the conditioning-direction (``x``) kernel (any ``\\int K = 1`` callable; default
  `gaussian_kernel`); as noted above it cancels in the ratio, so a symmetric kernel is almost always the
  right choice.
- `kernely` — the response-direction (``y``) kernel (default `gaussian_kernel`). A **symmetric** kernel
  gives ``\\int\\hat f\\,\\mathrm dy = 1`` exactly; an **asymmetric** `gamma_kernel` / `beta_kernel`, for a
  bounded / positive response, is boundary-bias-free but needs `normalize = true` for an exact integral —
  or supply your own [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel).
- `normalize` — when `true`, divide by ``\\int\\hat f(y\\mid x)\\,\\mathrm dy`` over the response range so
  the curve integrates to 1 **exactly**: the finishing step for an asymmetric `kernely`, a no-op for a
  symmetric one.

```jldoctest
julia> g = Histogram2D(-1.0, 1.0, 2, -1.0, 1.0, 2);

julia> observe_pairs!(g, [0.4], [0.4]);   # one pair → cell with centers (0.5, 0.5)

julia> round(conditional_density(g, 0.5, 0.5; hx = 0.5, hy = 0.2); digits = 4)
1.9947
```

# References

- M. Rosenblatt (1969), [*Conditional probability density and regression estimators*](https://openlibrary.org/works/OL7299345W) (no DOI), in *Multivariate Analysis II*, Academic Press.
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

**Estimand.** The **two-time transition density**
``p_\\tau(y \\mid x) = f_{X_{t+\\tau} \\mid X_t}(y \\mid x)`` — the conditional density of the state a
lag ``\\tau`` later, given the present state (nonstationary: ``p_{t,\\tau}``). For the
Ornstein–Uhlenbeck process this is Gaussian in closed form.

**Estimator.** The [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)
double-kernel estimator applied to the per-lag joint histogram of pairs
``\\bigl(X_t,  X_{t+\\tau}\\bigr)`` (binned by
[`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!) at the summary's lag
``\\tau``); the computation is identical — the name marks the two-time reading. A vector `y` returns the
transition-density curve at the fixed present state `x`.

# Arguments

- `summary` — the per-lag joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X_t, X_{t+\\tau})``.
- `x` — the present state, the conditioning value (a scalar).
- `y` — the later state(s) at which the density is evaluated (a scalar, or a vector for the curve).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``), as for
  `conditional_density`.
- `kernelx`, `kernely` — the conditioning and response kernels, as for `conditional_density` (an asymmetric
  `kernely` for a bounded / positive state).
- `normalize` — as for `conditional_density`: `true` makes ``\\int \\hat p_\\tau(y\\mid x)\\,\\mathrm dy = 1``
  exact for an asymmetric `kernely`.
"""
transition_density(summary::Histogram2D, x::Real, y; hx::Real, hy::Real,
                   kernelx = gaussian_kernel, kernely = gaussian_kernel, normalize::Bool = false) =
  conditional_density(summary, x, y; hx = hx, hy = hy, kernelx = kernelx, kernely = kernely,
                      normalize = normalize)


# The double-kernel variance R(K₁)R(K₂)f/(n h₁ h₂ f_X) is for symmetric (translation) kernels; the
# asymmetric gamma/beta kernels have a different variance structure — a recorded gap, as for the marginal
# density SE.
function _check_cd_se(kernelx, kernely)
  (kernelx isa AsymmetricKernel || kernely isa AsymmetricKernel) && throw(ArgumentError(
    "the conditional-density standard error is for symmetric kernels; the asymmetric gamma/beta-kernel variance is not supported"))
end


# The conditional-density estimate and its pointwise standard error at (x, y) from the pre-reduced X-profile
# (P, D): the two-term double-kernel variance R(K₁)R(K₂)·f̂/(n h₁ h₂ f̂_X(x)) − f̂²/n, with n f̂_X(x) = D the
# X-window denominator `_y_profile` already forms (Hyndman, Bashtannyk & Grunwald 1996; the −f̂²/n finite-
# sample term is the double-kernel analog of the marginal Wand–Jones correction), so
# SE = √(R(K₁)R(K₂)·f̂/(h₁ h₂ D) − f̂²/n). Both NaN where no data is local to x (D = 0). Shared by the SE
# reader and the paired value-plus-error reader, so the two agree by construction.
function _cd_value_and_se(summary::Histogram2D, P, D, y::Real, hx::Real, hy::Real, kernelx, kernely)
  f = _density_at(summary, P, D, y, hy, kernely)
  isfinite(f) || return (f, oftype(f, NaN))
  roughness = _quad(u -> kernelx(u)^2) * _quad(u -> kernely(u)^2)   # R(K₁)·R(K₂)
  fnn = max(f, 0.0)
  variance = roughness * fnn / (hx * hy * D) - fnn^2 / _total(summary)
  return (f, sqrt(max(variance, 0.0)))
end


"""
    conditional_density_standard_error(summary::Histogram2D, x, y;
                                       hx, hy, kernelx = gaussian_kernel,
                                       kernely = gaussian_kernel)   # ≡ SE[f̂(y∣x)]

**Estimand.** The pointwise standard error of the conditional-density estimator
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density),
``\\mathrm{SE}(y\\mid x) = \\sqrt{\\operatorname{Var}\\hat f(y\\mid x)}`` with the double-kernel asymptotic
variance

```math
\\operatorname{Var}\\hat f(y\\mid x) \\approx
  \\frac{R(K_1) R(K_2) f_{Y\\mid X}(y\\mid x)}{n h_1 h_2 f_X(x)}
  - \\frac{f_{Y\\mid X}(y\\mid x)^2}{n},
  \\qquad R(K) = \\int K^2 .
```

Here ``f_X`` is the marginal density of the conditioning variable ``X``, and the ``-f^2/n`` term is the
double-kernel analog of the marginal Wand–Jones correction.

**Estimator.** The two-term ``\\sqrt{R(K_1)R(K_2) \\hat f(y\\mid x)/(h_1 h_2 D) - \\hat f(y\\mid x)^2/n}``
from the recovered conditional
density ``\\hat f(y\\mid x)``, the bandwidths ``h_1`` and ``h_2`` (`hx` and `hy`), the kernel
roughnesses ``R(K_1)``, ``R(K_2)`` (by quadrature), and the **``x``-window denominator**
``D = \\sum_j N_j K_{h_1}(c^x_j - x) = n \\hat f_X(x)`` — the very sum `conditional_density` forms to
normalize its profile, so ``n \\hat f_X(x)`` costs nothing extra. It is the two-dimensional analog of the
marginal [`density_standard_error`](@ref GBGStreamingStatistics.MarginalDensity.density_standard_error),
with the ``x``-window effective sample size ``n h_1 \\hat f_X(x)`` in place of ``n``; the ``z``-band
``\\hat f(y\\mid x) \\pm z \\mathrm{SE}`` covers ``f_{Y\\mid X}`` where the ``O(h^2)`` smoothing bias is
small. Symmetric kernels only (the asymmetric-`kernely` variance is a recorded gap). `NaN` where no data is
local to ``x``.

# Arguments

The arguments are as for `conditional_density`, without `normalize` (the band is for the un-normalized
estimate):

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X, Y)``.
- `x` — the conditioning value (a scalar).
- `y` — the response point(s) at which the standard error is evaluated (a scalar, or a vector).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``).
- `kernelx`, `kernely` — the conditioning and response kernels (default [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)); symmetric
  kernels only.

Reference: R. J. Hyndman, D. M. Bashtannyk & G. K. Grunwald (1996), [*Estimating and visualizing
conditional densities*](https://doi.org/10.1080/10618600.1996.10474715), J. Comput. Graph. Statist. 5(4).
"""
function conditional_density_standard_error(summary::Histogram2D, x::Real, y::Real; hx::Real, hy::Real,
                                            kernelx = gaussian_kernel, kernely = gaussian_kernel)
  _check_bandwidths(hx, hy)
  _check_cd_se(kernelx, kernely)
  P, D = _y_profile(summary, x, hx, kernelx)
  return _cd_value_and_se(summary, P, D, y, hx, hy, kernelx, kernely)[2]
end

function conditional_density_standard_error(summary::Histogram2D, x::Real, ys::AbstractVector; hx::Real,
                                            hy::Real, kernelx = gaussian_kernel, kernely = gaussian_kernel)
  _check_bandwidths(hx, hy)
  _check_cd_se(kernelx, kernely)
  P, D = _y_profile(summary, x, hx, kernelx)
  return [_cd_value_and_se(summary, P, D, y, hx, hy, kernelx, kernely)[2] for y in ys]
end


"""
    conditional_density_with_standard_error(summary::Histogram2D, x, y;
                                            hx, hy, kernelx = gaussian_kernel,
                                            kernely = gaussian_kernel)

The estimate and its standard error from one ``x``-profile reduction — the `NamedTuple`
``(\\texttt{value}, \\texttt{standard\\_error})``, with `value` the
[`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)
``\\hat f(y\\mid x)`` and `standard_error` its
[`conditional_density_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_standard_error).
They share the one profile, so they describe the same estimator by construction. For a vector `y`, returns
a `NamedTuple` of vectors.

# Arguments

The arguments are as for `conditional_density`, without `normalize` (the band is for the un-normalized
estimate):

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X, Y)``.
- `x` — the conditioning value (a scalar).
- `y` — the response point(s) (a scalar, or a vector for the curve and its band).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``).
- `kernelx`, `kernely` — the conditioning and response kernels (default [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)); symmetric
  kernels only.
"""
function conditional_density_with_standard_error(summary::Histogram2D, x::Real, y::Real; hx::Real, hy::Real,
                                                 kernelx = gaussian_kernel, kernely = gaussian_kernel)
  _check_bandwidths(hx, hy)
  _check_cd_se(kernelx, kernely)
  P, D = _y_profile(summary, x, hx, kernelx)
  v, se = _cd_value_and_se(summary, P, D, y, hx, hy, kernelx, kernely)
  return (value = v, standard_error = se)
end

function conditional_density_with_standard_error(summary::Histogram2D, x::Real, ys::AbstractVector; hx::Real,
                                                 hy::Real, kernelx = gaussian_kernel, kernely = gaussian_kernel)
  _check_bandwidths(hx, hy)
  _check_cd_se(kernelx, kernely)
  P, D = _y_profile(summary, x, hx, kernelx)
  vses = [_cd_value_and_se(summary, P, D, y, hx, hy, kernelx, kernely) for y in ys]
  return (value = [vs[1] for vs in vses], standard_error = [vs[2] for vs in vses])
end


"""
    transition_density_standard_error(summary::Histogram2D, x, y;
                                      hx, hy, kernelx = gaussian_kernel,
                                      kernely = gaussian_kernel)

The standard error of [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density)
— identical to [`conditional_density_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_standard_error)
on the per-lag joint histogram (the name marks the two-time reading), paired with the estimate by
[`transition_density_with_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density_with_standard_error).

# Arguments

The arguments are as for `transition_density`, without `normalize` (the band is for the un-normalized
estimate):

- `summary` — the per-lag joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X_t, X_{t+\\tau})``.
- `x` — the present state, the conditioning value (a scalar).
- `y` — the later state(s) at which the standard error is evaluated (a scalar, or a vector).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``).
- `kernelx`, `kernely` — the conditioning and response kernels (default [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)); symmetric
  kernels only.
"""
transition_density_standard_error(summary::Histogram2D, x::Real, y; hx::Real, hy::Real,
                                  kernelx = gaussian_kernel, kernely = gaussian_kernel) =
  conditional_density_standard_error(summary, x, y; hx = hx, hy = hy, kernelx = kernelx, kernely = kernely)


"""
    transition_density_with_standard_error(summary::Histogram2D, x, y;
                                           hx, hy, kernelx = gaussian_kernel,
                                           kernely = gaussian_kernel)

The transition density
[`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density) and its standard
error from one pass (the `(value, standard_error)` `NamedTuple`) — the two-time reading of
[`conditional_density_with_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_with_standard_error).

# Arguments

The arguments are as for `transition_density`, without `normalize` (the band is for the un-normalized
estimate):

- `summary` — the per-lag joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X_t, X_{t+\\tau})``.
- `x` — the present state, the conditioning value (a scalar).
- `y` — the later state(s) (a scalar, or a vector for the curve and its band).
- `hx`, `hy` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``).
- `kernelx`, `kernely` — the conditioning and response kernels (default [`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel)); symmetric
  kernels only.
"""
transition_density_with_standard_error(summary::Histogram2D, x::Real, y; hx::Real, hy::Real,
                                       kernelx = gaussian_kernel, kernely = gaussian_kernel) =
  conditional_density_with_standard_error(summary, x, y; hx = hx, hy = hy, kernelx = kernelx, kernely = kernely)

end # module ConditionalDensity
