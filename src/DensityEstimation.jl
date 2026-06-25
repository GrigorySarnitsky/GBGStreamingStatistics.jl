# DensityEstimation — §4.1: the marginal density f_X, estimated from a 1-D histogram by binned KDE
# (WARPing). Written estimand-first (`docs/ARCHITECTURE.md` §2, §4.1).

module DensityEstimation

using ..Histograms: Histogram1D, _nbins, _bincenter, _total, _sample_std

export density, adaptive_density, gaussian_kernel, epanechnikov_kernel


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


# Silverman pilot bandwidth from the recovered marginal ``\hat\sigma``.
_silverman(summary::Histogram1D) =
  (4 / 3)^(1 / 5) * _sample_std(summary) * float(_total(summary))^(-1 / 5)


# Abramson per-bin bandwidths: a fixed-bandwidth pilot ``\tilde f`` at h0 gives ``h_j = h0 (\tilde
# f(c_j)/g)^{-1/2}`` with ``g`` the geometric mean of ``\tilde f`` over the sample (occupied bins,
# weighted by counts). Host-side — the O(nbins²) pilot brings the counts to the host.
function _adaptive_bandwidths(summary::Histogram1D, h0, kernel)
  n = _total(summary)
  nb = _nbins(summary)
  counts = Array(summary.counts)
  centers = [_bincenter(summary, j) for j in 1:nb]
  invh0 = inv(h0)
  pilot = zeros(Float64, nb)
  @inbounds for j in 1:nb
    iszero(counts[j]) && continue
    cj = centers[j]
    acc = 0.0
    for k in 1:nb
      ck = counts[k]
      iszero(ck) && continue
      acc += ck * kernel((cj - centers[k]) * invh0)
    end
    pilot[j] = acc / (n * h0)
  end
  logg = 0.0
  @inbounds for j in 1:nb
    (counts[j] > 0 && pilot[j] > 0) || continue
    logg += counts[j] * log(pilot[j])
  end
  g = exp(logg / n)
  hj = [pilot[j] > 0 ? h0 * sqrt(g / pilot[j]) : h0 for j in 1:nb]
  return counts, centers, hj
end


function _adaptive_eval(counts, centers, hj, n, x, kernel)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel((x - centers[j]) / hj[j]) / hj[j]
  end
  return acc / n
end


"""
    adaptive_density(summary::Histogram1D, x; h0 = nothing, kernel = gaussian_kernel)   # ≡ f_X

**Estimand.** the marginal density ``f_X`` — the same population quantity as
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density).

**Estimator** (Abramson variable-bandwidth / sample-point KDE). A fixed-bandwidth pilot ``\\tilde f``
at ``h_0`` sets a **per-bin bandwidth** ``h_j = h_0\\,(\\tilde f(c_j)/g)^{-1/2}`` (``g`` the geometric
mean of ``\\tilde f`` over the sample), then

```math
\\hat f(x) = \\frac1n \\sum_j \\frac{n_j}{h_j}\\, K\\!\\left(\\frac{x - c_j}{h_j}\\right).
```

The kernel **widens where the data is sparse** (a small pilot ``\\tilde f(c_j)`` gives a large
``h_j``) and narrows in the peak, so the **tails** carry far less variance than a fixed bandwidth — the
regime that matters for heavy / intermittent distributions. Abramson's square-root law lowers the bias
to ``O(h^4)`` (from the ``O(h^2)`` of
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density)). The pilot ``h_0`` defaults to the
global Silverman rule; `kernel` is any ``K`` with ``\\int K = 1``.

Reference: I. S. Abramson (1982), [*On bandwidth variation in kernel estimates—a square root
law*](https://doi.org/10.1214/aos/1176345986), Ann. Statist. 10(4).
"""
function adaptive_density(summary::Histogram1D, x::Real; h0 = nothing, kernel = gaussian_kernel)
  n = _total(summary)
  iszero(n) && return 0.0
  h = h0 === nothing ? _silverman(summary) : float(h0)
  h > 0 || throw(ArgumentError("adaptive_density needs a positive pilot bandwidth (got h0 = $h0; the Silverman default vanishes with < 2 spread-out samples — supply h0)"))
  counts, centers, hj = _adaptive_bandwidths(summary, h, kernel)
  return _adaptive_eval(counts, centers, hj, n, x, kernel)
end

function adaptive_density(summary::Histogram1D, xs::AbstractVector; h0 = nothing, kernel = gaussian_kernel)
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  h = h0 === nothing ? _silverman(summary) : float(h0)
  h > 0 || throw(ArgumentError("adaptive_density needs a positive pilot bandwidth (got h0 = $h0)"))
  counts, centers, hj = _adaptive_bandwidths(summary, h, kernel)
  return [_adaptive_eval(counts, centers, hj, n, x, kernel) for x in xs]
end

end # module DensityEstimation
