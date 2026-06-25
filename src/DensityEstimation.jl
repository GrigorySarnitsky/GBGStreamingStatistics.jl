# DensityEstimation — §4.1: the marginal density f_X, estimated from a 1-D histogram by binned KDE
# (WARPing). Written estimand-first (`docs/ARCHITECTURE.md` §2, §4.1).

module DensityEstimation

using ..Histograms: Histogram1D, _nbins, _bincenter, _total

export density, gaussian_kernel, epanechnikov_kernel


"""
    gaussian_kernel(u)

The standard Gaussian smoothing kernel ``K(u) = e^{-u^2/2}/\\sqrt{2\\pi}`` (``\\int K = 1``,
``\\mu_2(K) = 1``, roughness ``R(K) = 1/(2\\sqrt{\\pi})``) — the
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density) default.
"""
@inline gaussian_kernel(u) = exp(-u^2 / 2) / sqrt(2π)


"""
    epanechnikov_kernel(u)

The **Epanechnikov** kernel ``K(u) = \\tfrac{3}{4}(1 - u^2)`` for ``|u| \\le 1`` (zero outside) — the
MSE-optimal second-order kernel (``\\int K = 1``, ``\\mu_2(K) = 1/5``, roughness ``R(K) = 3/5``). Pass
it as `kernel = epanechnikov_kernel`; its compact support touches only the bins within ``h`` of the
query point. Note [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) returns a
**Gaussian**-tuned ``h`` — the AMISE-optimal Epanechnikov bandwidth is ``\\approx 2.21\\times`` larger,
so scale it (or supply ``h`` directly).
"""
@inline epanechnikov_kernel(u) = abs(u) < 1 ? 0.75 * (1 - u^2) : 0.0


"""
    density(summary::Histogram1D, x; h, kernel = gaussian_kernel)   # ≡ f_X

**Estimand.** the marginal (Lebesgue) density ``f_X = \\mathrm{d}\\mu_X/\\mathrm{d}\\lambda`` of the
conditioning variable ``X``, defined by ``P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda``.

**Estimator** (binned KDE / WARPing). From the bin counts ``n_j`` at centers ``c_j`` of `summary`
(a [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D)), the total count ``n``, and a
kernel ``K`` with ``\\int K = 1`` and bandwidth ``h``,

```math
\\hat f_h(x) = \\frac{1}{n\\,h} \\sum_{j} n_j\\, K\\!\\left(\\frac{x - c_j}{h}\\right),
```

the binned approximation to ``\\frac{1}{nh}\\sum_i K\\!\\left(\\frac{x - X_i}{h}\\right)`` ([Hall &
Wand 1996](https://doi.org/10.1006/jmva.1996.0009)). The bandwidth ``h`` comes from
[`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth) (or is supplied directly); `kernel`
defaults to [`gaussian_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gaussian_kernel) — pass
[`epanechnikov_kernel`](@ref GBGStreamingStatistics.DensityEstimation.epanechnikov_kernel), or any ``K``
with ``\\int K = 1``. Here the convolution is a direct sum over occupied bins; the
FFT amortization of the sum over a whole bandwidth grid is forthcoming (`CLAUDE.md` roadmap). A
vector of query points returns a vector of densities.

```jldoctest
julia> h = Histogram1D(-1.0, 1.0, 2);

julia> observe!(h, [0.4]);          # one sample, binned at center 0.5

julia> round(density(h, 0.5; h = 0.2); digits = 4)
1.9947
```
"""
function density(summary::Histogram1D, x::Real; h::Real, kernel = gaussian_kernel)
  h > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  n = _total(summary)
  iszero(n) && return 0.0
  invh = inv(float(h))
  acc = 0.0
  @inbounds for j in 1:_nbins(summary)
    c = summary.counts[j]
    iszero(c) && continue
    acc += c * kernel((x - _bincenter(summary, j)) * invh)
  end
  return acc / (n * h)
end

density(summary::Histogram1D, xs::AbstractVector; kwargs...) =
  [density(summary, x; kwargs...) for x in xs]

end # module DensityEstimation
