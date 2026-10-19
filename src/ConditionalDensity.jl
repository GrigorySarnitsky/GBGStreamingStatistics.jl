# ConditionalDensity — §4.3: the conditional density f_{Y|X} and the two-time transition density
# p_τ(y|x), from a joint histogram by the double-kernel estimator (`docs/ARCHITECTURE.md` §4.3).

module ConditionalDensity

using ..Histograms: Histogram2D, _xbincenter, _ybincenter, _total, _host
using ..Kernels: gaussian_kernel, _kernel_weight, AsymmetricKernel, _quad

export conditional_density, conditional_density_standard_error, conditional_density_with_standard_error,
       transition_density, transition_density_standard_error, transition_density_with_standard_error


_check_bandwidths(h_x::Real, h_y::Real) =
  (h_x > 0 && h_y > 0) ||
    throw(ArgumentError("bandwidths h_x, h_y must be positive (got h_x = $h_x, h_y = $h_y)"))


# Reduce the joint histogram at conditioning value `x` to a 1-D x-kernel-weighted y-profile:
# P_l = Σ_j n_{jl} K_{h_x}(c^x_j, x) (the K_{h_x}-weighted count in y-row l) and D = Σ_l P_l = Σ_j N_j
# K_{h_x}(c^x_j, x) (the marginal-X KDE numerator). The conditional distribution over y-cells is P/D;
# smoothing it in y gives the density. The per-column weight is the shared `_kernel_weight` atom (symmetric
# or asymmetric `kernel_x` alike); the x-kernel's normalization cancels in the f_{X,Y}/f_X ratio, so it is
# irrelevant here. Host-side (`_host`: no copy for a CPU summary).
function _y_profile(summary::Histogram2D, x::Real, kernel_x::KX, h_x::Real) where {KX}
  counts = _host(summary.counts)
  n_x = summary.n_x
  n_y = summary.n_y
  gx = Vector{Float64}(undef, n_x)
  @inbounds for j in 1:n_x
    gx[j] = _kernel_weight(_xbincenter(summary, j), x, kernel_x, h_x)
  end
  P = zeros(Float64, n_y)
  @inbounds for k in 1:n_y
    base = (k - 1) * n_x
    acc = 0.0
    for j in 1:n_x
      c = counts[base + j]
      iszero(c) && continue
      acc += c * gx[j]
    end
    P[k] = acc
  end
  return P, sum(P)
end


# Smooth the y-profile to the density at `y`. The double-kernel ratio f_{X,Y}/f_X cancels the x-kernel's
# normalization (it appears in both the profile P and its total D); only the y-direction normalization
# survives, carried by the shared `_kernel_weight` atom (its 1/h_y folded in). For a symmetric `kernel_y`,
# ∫K = 1 makes ∫ f̂_{Y|X}(y | x) dy = 1 exactly; an asymmetric `kernel_y` (a normalized density in its data
# argument) gives ∫ f̂_{Y|X} dy ≈ 1 — exact only after `_response_normalizer` (as for the marginal
# `density`).
function _density_at(summary::Histogram2D, P, D, y::Real, kernel_y::KY, h_y::Real) where {KY}
  D > 0 || return oftype(D, NaN)
  acc = 0.0
  @inbounds for k in 1:summary.n_y
    Pk = P[k]
    iszero(Pk) && continue
    acc += Pk * _kernel_weight(_ybincenter(summary, k), y, kernel_y, h_y)
  end
  return acc / D
end


# ∫ f̂_{Y|X}(y | x) dy over the summary's y-range ∩ the response kernel's support, by the trapezoidal rule
# on a fine grid — the per-x renormalization that makes ∫ f̂_{Y|X} dy = 1 exactly for a non-self-normalizing
# asymmetric response kernel (a symmetric kernel already gives 1, so `normalize` is a no-op there).
function _response_normalizer(summary::Histogram2D, P, D, kernel_y::KY, h_y::Real) where {KY}
  lower, upper = float(summary.lower_y), float(summary.upper_y)
  if kernel_y isa AsymmetricKernel
    lower = max(lower, kernel_y.support[1])
    upper = min(upper, kernel_y.support[2])
  end
  npts = 1001
  grid = range(lower, upper; length = npts)
  acc = 0.0
  @inbounds for (i, yy) in enumerate(grid)
    w = (i == 1 || i == npts) ? 0.5 : 1.0
    acc += w * _density_at(summary, P, D, yy, kernel_y, h_y)
  end
  return acc * step(grid)
end


"""
    conditional_density(summary::Histogram2D, x, y; kernel_x = gaussian_kernel,
                        kernel_y = gaussian_kernel, h_x, h_y, normalize = false)

**Estimand.** The conditional density
``f_{Y\\mid X}(y \\mid x) = f_{X,Y}(x, y) / f_X(x)``.

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996). From the joint
[`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) `summary`, with a conditioning kernel
``K_x`` (`kernel_x`, of bandwidth ``h_x`` given by `h_x`, in ``x``) and a response kernel ``K_y`` (`kernel_y`,
of bandwidth ``h_y`` given by `h_y`, in ``y``); both satisfy ``\\int K = 1`` and default to
[`gaussian_kernel`](@ref GBGStreamingStatistics.Kernels.gaussian_kernel). Over the observed pairs of
conditioning value and response, as many as the sample size, the estimator is

```math
\\hat f_{Y\\mid X}(y \\mid x) =
  \\frac{\\sum_i K_{h_x}(X_i - x)  K_{h_y}(Y_i - y)}{\\sum_i K_{h_x}(X_i - x)} ,
```

evaluated from the bin counts ``n_{jl}`` at cell centers ``(c^x_j, c^y_l)`` (the binned form, with
``N_j = \\sum_l n_{jl}`` the ``x``-column count):

```math
\\hat f_{Y\\mid X}(y \\mid x) =
  \\frac{\\sum_{j,l} n_{jl}  K_{h_x}(c^x_j - x)  K_{h_y}(c^y_l - y)}
       {\\sum_{j} N_j  K_{h_x}(c^x_j - x)} .
```

Only the ``x``-kernel's normalization cancels in the ``f_{X,Y}/f_X`` ratio; its weighting sets the conditioning window, a degree-0 (Nadaraya–Watson-type)
weighting that keeps a first-order bias at a bounded conditioning edge.

For a symmetric response `kernel_y` the surviving ``1/h_y`` makes ``\\int \\hat f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y = 1``
exactly. For a bounded / positive response pass an asymmetric `kernel_y`
([`gamma_kernel`](@ref GBGStreamingStatistics.Kernels.gamma_kernel) for ``[0,\\infty)``,
[`beta_kernel`](@ref GBGStreamingStatistics.Kernels.beta_kernel) for ``[0,1]``), the boundary-bias fix the
marginal [`density`](@ref GBGStreamingStatistics.MarginalDensity.density) uses: a density
``K_y(c^y_l; y, b)`` in its data argument (the cell center), of scale ``b`` set by `h_y`, so no mass
leaks past the edge, but only
approximately normalized in ``y``; set `normalize = true` for ``\\int \\hat f_{Y\\mid X}\\,\\mathrm{d}y = 1`` exactly.

The conditional MISE falls as ``O(n^{-2/3})`` in the sample size, slower than a marginal density's. A vector `y` returns the
curve at the fixed conditioning value `x`; where no data is local to `x` the
estimate is `NaN`.

# Arguments

- `summary` — the joint `Histogram2D` of ``(X, Y)``.
- `x` — the conditioning value (a scalar).
- `y` — the response point(s) at which the density is evaluated (scalar, or a vector for the curve).
- `kernel_x` — the ``x``-direction kernel (default `gaussian_kernel`); its normalization cancels in the
  ratio; its weighting sets the window.
- `kernel_y` — the ``y``-direction kernel (default `gaussian_kernel`): symmetric, with
  ``\\int\\hat f_{Y\\mid X}\\,\\mathrm dy = 1`` exactly, or an asymmetric `gamma_kernel` / `beta_kernel` (or your
  own [`AsymmetricKernel`](@ref GBGStreamingStatistics.Kernels.AsymmetricKernel)) for a bounded / positive
  response, boundary-bias-free but needing `normalize = true`.
- `h_x`, `h_y` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``), supplied.
- `normalize` — when `true`, divide by ``\\int\\hat f_{Y\\mid X}(y\\mid x)\\,\\mathrm dy`` over the response range so
  the curve integrates to 1 exactly: needed for an asymmetric `kernel_y`, a no-op for a symmetric one.

```jldoctest
julia> joint = Histogram2D(-1.0, 1.0, 2, -1.0, 1.0, 2);

julia> # one pair → cell with centers (0.5, 0.5)
       observe_pairs!(joint, [0.4], [0.4]);

julia> round(conditional_density(joint, 0.5, 0.5; h_x = 0.5, h_y = 0.2); digits = 4)
1.9947
```

# References

- M. Rosenblatt (1969), [*Conditional probability density and regression estimators*](https://openlibrary.org/works/OL7299345W) (no DOI), in *Multivariate Analysis II*, Academic Press.
- R. J. Hyndman, D. M. Bashtannyk & G. K. Grunwald (1996), [*Estimating and visualizing conditional densities*](https://doi.org/10.1080/10618600.1996.10474715), J. Comput. Graph. Statist. 5(4).
"""
function conditional_density(summary::Histogram2D, x::Real, y::Real;
                             kernel_x = gaussian_kernel, kernel_y = gaussian_kernel,
                             h_x::Real, h_y::Real, normalize::Bool = false)
  _check_bandwidths(h_x, h_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  f = _density_at(summary, P, D, y, kernel_y, h_y)
  (normalize && D > 0) || return f
  return f / _response_normalizer(summary, P, D, kernel_y, h_y)
end

function conditional_density(summary::Histogram2D, x::Real, ys::AbstractVector;
                             kernel_x = gaussian_kernel, kernel_y = gaussian_kernel,
                             h_x::Real, h_y::Real, normalize::Bool = false)
  _check_bandwidths(h_x, h_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  Z = (normalize && D > 0) ? _response_normalizer(summary, P, D, kernel_y, h_y) : 1.0
  return [_density_at(summary, P, D, y, kernel_y, h_y) / Z for y in ys]
end


"""
    transition_density(summary::Histogram2D, x, y; kernel_x = gaussian_kernel,
                       kernel_y = gaussian_kernel, h_x, h_y, normalize = false)

**Estimand.** The two-time transition density
``p_\\tau(y \\mid x) = f_{Y(t+\\tau) \\mid X(t)}(y \\mid x)`` — the conditional density of the state a
lag ``\\tau`` later, given the present state (``p_{t,\\tau}`` if it also depends on
time). For the
Ornstein–Uhlenbeck process this is Gaussian in closed form.

**Estimator.** The [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density)
double-kernel estimator applied to a joint histogram fed the pairs ``\\bigl(X(t), Y(t+\\tau)\\bigr)`` at one
fixed lag ``\\tau`` (binned by
[`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!); the summary does not record
the lag); the computation is identical — the name marks the two-time reading. A vector `y` returns the
transition-density curve at the fixed present state `x`.

# Arguments

- `summary` — the per-lag joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of ``(X(t), Y(t+\\tau))``.
- `x` — the present state, the conditioning value (a scalar).
- `y` — the later state(s) at which the density is evaluated (a scalar, or a vector for the curve).
- `kernel_x`, `kernel_y` — the conditioning and response kernels, as for `conditional_density` (an asymmetric
  `kernel_y` for a bounded / positive state).
- `h_x`, `h_y` — the conditioning (``x``) and response (``y``) bandwidths (both ``> 0``), as for
  `conditional_density`.
- `normalize` — as for `conditional_density`: `true` makes ``\\int \\hat p_\\tau(y\\mid x)\\,\\mathrm dy = 1``
  exact for an asymmetric `kernel_y`.
"""
transition_density(summary::Histogram2D, x::Real, y; kernel_x = gaussian_kernel,
                   kernel_y = gaussian_kernel, h_x::Real, h_y::Real,
                   normalize::Bool = false) =
  conditional_density(summary, x, y; h_x = h_x, h_y = h_y, kernel_x = kernel_x, kernel_y = kernel_y,
                      normalize = normalize)


# The double-kernel variance R(K_x)R(K_y) f_{Y|X}/(n h_x h_y f_X) is for symmetric (translation) kernels;
# the asymmetric gamma/beta kernels have a different variance structure — a recorded gap, as for the
# marginal density SE.
function _check_cd_se(kernel_x, kernel_y)
  (kernel_x isa AsymmetricKernel || kernel_y isa AsymmetricKernel) && throw(ArgumentError(
    "the conditional-density standard error is for symmetric kernels; the asymmetric gamma/beta-kernel variance is not supported"))
end


# The conditional-density estimate and its pointwise standard error at (x, y) from the pre-reduced
# x-profile (P, D). The estimate is a weighted mean of K_{h_y}(Y − y) over the samples near x, so its variance
# is the sum of the squared weights, R(K_x)/(h_x D), times the variance of one response-kernel term,
# R(K_y) f̂_{Y|X}/h_y − f̂_{Y|X}² (Hyndman, Bashtannyk & Grunwald 1996), with n f̂_X(x) = D the x-window
# denominator `_y_profile` already forms: SE = √(R(K_x)/(h_x D) · (R(K_y) f̂_{Y|X}/h_y − f̂_{Y|X}²)). The
# second term carries the same local weight as the first — the marginal −f̂_X²/n has weight 1/n because
# every sample counts there, which is not so here. Both NaN where no data is local to x (D = 0). Shared by
# the SE reader and the paired value-plus-error reader, so the two agree by construction.
function _cd_value_and_se(summary::Histogram2D, P, D, y::Real, kernel_x::KX, kernel_y::KY, h_x::Real, h_y::Real) where {KX, KY}
  f = _density_at(summary, P, D, y, kernel_y, h_y)
  isfinite(f) || return (f, oftype(f, NaN))
  roughness_x = _quad(u -> kernel_x(u)^2)                       # R(K_x)
  roughness_y = _quad(u -> kernel_y(u)^2)                       # R(K_y)
  fnn = max(f, 0.0)
  variance = roughness_x / (h_x * D) * (roughness_y * fnn / h_y - fnn^2)
  return (f, sqrt(max(variance, 0.0)))
end


"""
    conditional_density_standard_error(summary::Histogram2D, x, y;
                                       kernel_x = gaussian_kernel,
                                       kernel_y = gaussian_kernel, h_x, h_y)

**Estimand.** The pointwise standard error of the conditional-density estimator
`conditional_density`,
``\\operatorname{SE}(y\\mid x) = \\sqrt{\\operatorname{Var}\\hat f_{Y\\mid X}(y\\mid x)}`` with the double-kernel asymptotic
variance

```math
\\operatorname{Var}\\hat f_{Y\\mid X}(y\\mid x) \\approx
  \\frac{R(K_x)}{n h_x f_X(x)}
  \\biggl[\\frac{R(K_y) f_{Y\\mid X}(y\\mid x)}{h_y} - f_{Y\\mid X}(y\\mid x)^2\\biggr],
  \\qquad R(K) = \\int K^2 .
```

Here ``f_X`` is the marginal density of ``X``; the bracket is the variance of one response-kernel term
``K_{h_y}(Y - y)``, and ``R(K_x)/(n h_x f_X(x))`` the sum of the squared weights of the samples near ``x``.

**Estimator.** It is the two-term ``\\sqrt{R(K_x)/(h_x D) \\cdot [R(K_y) \\hat f_{Y\\mid X}(y\\mid x)/h_y - \\hat f_{Y\\mid X}(y\\mid x)^2]}``
from the conditional density
``\\hat f_{Y\\mid X}(y\\mid x)``, the bandwidths ``h_x`` (`h_x`) and ``h_y`` (`h_y`), the kernel
roughnesses ``R(K_x)``, ``R(K_y)``, and the ``x``-window denominator
``D = \\sum_j N_j K_{h_x}(c^x_j - x) = n \\hat f_X(x)``, the sum `conditional_density`
normalizes by. It is the two-dimensional analog of
[`density_standard_error`](@ref GBGStreamingStatistics.MarginalDensity.density_standard_error), with the ``x``-window effective sample size ``n h_x \\hat f_X(x)`` in place of the sample size ``n``; the 95% band
``\\hat f_{Y\\mid X}(y\\mid x) \\pm 1.96 \\operatorname{SE}`` covers ``f_{Y\\mid X}`` where the ``O(h^2)`` bias, in
either bandwidth, is small. It takes symmetric kernels only (no asymmetric-`kernel_y` variance yet) and is `NaN` where no data is local to ``x``.

Arguments: as for [`conditional_density`](@ref
GBGStreamingStatistics.ConditionalDensity.conditional_density), without `normalize` (the band is for the
un-normalized estimate).

Reference: R. J. Hyndman, D. M. Bashtannyk & G. K. Grunwald (1996), [*Estimating and visualizing
conditional densities*](https://doi.org/10.1080/10618600.1996.10474715), J. Comput. Graph. Statist. 5(4).
"""
function conditional_density_standard_error(summary::Histogram2D, x::Real, y::Real;
                                            kernel_x = gaussian_kernel,
                                            kernel_y = gaussian_kernel, h_x::Real,
                                            h_y::Real)
  _check_bandwidths(h_x, h_y)
  _check_cd_se(kernel_x, kernel_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  return _cd_value_and_se(summary, P, D, y, kernel_x, kernel_y, h_x, h_y)[2]
end

function conditional_density_standard_error(summary::Histogram2D, x::Real, ys::AbstractVector;
                                            kernel_x = gaussian_kernel,
                                            kernel_y = gaussian_kernel, h_x::Real,
                                            h_y::Real)
  _check_bandwidths(h_x, h_y)
  _check_cd_se(kernel_x, kernel_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  return [_cd_value_and_se(summary, P, D, y, kernel_x, kernel_y, h_x, h_y)[2] for y in ys]
end


"""
    conditional_density_with_standard_error(summary::Histogram2D, x, y;
                                            kernel_x = gaussian_kernel,
                                            kernel_y = gaussian_kernel, h_x, h_y)

Returns the estimate and its standard error from one ``x``-profile reduction, the `NamedTuple`
``(\\texttt{value}, \\texttt{standard\\_error})``, with `value` the
`conditional_density`
``\\hat f_{Y\\mid X}(y\\mid x)`` and `standard_error` its
[`conditional_density_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_standard_error).
They share the one profile, so they describe the same estimator by construction. For a vector `y`, returns
a `NamedTuple` of vectors.

Arguments: as for [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density),
without `normalize` (the band is for the un-normalized estimate); symmetric kernels only.

"""
function conditional_density_with_standard_error(summary::Histogram2D, x::Real, y::Real;
                                                 kernel_x = gaussian_kernel,
                                                 kernel_y = gaussian_kernel, h_x::Real,
                                                 h_y::Real)
  _check_bandwidths(h_x, h_y)
  _check_cd_se(kernel_x, kernel_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  v, se = _cd_value_and_se(summary, P, D, y, kernel_x, kernel_y, h_x, h_y)
  return (value = v, standard_error = se)
end

function conditional_density_with_standard_error(summary::Histogram2D, x::Real, ys::AbstractVector;
                                                 kernel_x = gaussian_kernel,
                                                 kernel_y = gaussian_kernel, h_x::Real,
                                                 h_y::Real)
  _check_bandwidths(h_x, h_y)
  _check_cd_se(kernel_x, kernel_y)
  P, D = _y_profile(summary, x, kernel_x, h_x)
  vses = [_cd_value_and_se(summary, P, D, y, kernel_x, kernel_y, h_x, h_y) for y in ys]
  return (value = [vs[1] for vs in vses], standard_error = [vs[2] for vs in vses])
end


"""
    transition_density_standard_error(summary::Histogram2D, x, y;
                                      kernel_x = gaussian_kernel,
                                      kernel_y = gaussian_kernel, h_x, h_y)

The standard error of `transition_density` is
[`conditional_density_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_standard_error)
on the per-lag joint histogram (the name marks the two-time reading); it is paired with the estimate by
[`transition_density_with_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density_with_standard_error).

Arguments: as for [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density),
without `normalize` (the band is for the un-normalized estimate).

"""
transition_density_standard_error(summary::Histogram2D, x::Real, y;
                                  kernel_x = gaussian_kernel,
                                  kernel_y = gaussian_kernel, h_x::Real, h_y::Real) =
  conditional_density_standard_error(summary, x, y; h_x = h_x, h_y = h_y, kernel_x = kernel_x, kernel_y = kernel_y)


"""
    transition_density_with_standard_error(summary::Histogram2D, x, y;
                                           kernel_x = gaussian_kernel,
                                           kernel_y = gaussian_kernel, h_x, h_y)

Returns the transition density `transition_density` and its standard error from one pass (the
`(value, standard_error)` `NamedTuple`), the two-time reading of
[`conditional_density_with_standard_error`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density_with_standard_error).

Arguments: as for [`transition_density`](@ref GBGStreamingStatistics.ConditionalDensity.transition_density),
without `normalize` (the band is for the un-normalized estimate).

"""
transition_density_with_standard_error(summary::Histogram2D, x::Real, y;
                                       kernel_x = gaussian_kernel,
                                       kernel_y = gaussian_kernel, h_x::Real, h_y::Real) =
  conditional_density_with_standard_error(summary, x, y; h_x = h_x, h_y = h_y, kernel_x = kernel_x, kernel_y = kernel_y)

end # module ConditionalDensity
