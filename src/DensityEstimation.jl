# DensityEstimation — §4.1: the marginal density f_X. A single `density` spans the grid of a
# local-log-polynomial degree (0 = KDE, ≥1 = local likelihood) crossed with a fixed or Abramson-variable
# bandwidth. Two stages: set the bandwidth, then fit. The bandwidth is *sample-point* (per-bin h_j) for
# the KDE and *balloon* (one h(x) at the query) for the local-log fit. Written estimand-first
# (`docs/ARCHITECTURE.md` §2, §4.1).

module DensityEstimation

using LinearAlgebra

using ..Histograms: Histogram1D, _nbins, _bincenter, _total, _sample_std

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


# Silverman pilot/bandwidth from the recovered marginal ``\hat\sigma`` — the default `h`.
_silverman(summary::Histogram1D) =
  (4 / 3)^(1 / 5) * _sample_std(summary) * float(_total(summary))^(-1 / 5)


# Fixed-bandwidth binned KDE at x (single h) — the degree-0 fixed estimator and the adaptive pilot.
function _kde_h(counts, centers, n, x, h, kernel)
  invh = inv(h)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel((x - centers[j]) * invh)
  end
  return acc / (n * h)
end


# Abramson: per-bin (sample-point) bandwidths ``h_j = h (\tilde f(c_j)/g)^{-1/2}`` and the geometric
# mean ``g`` of the fixed-h pilot ``\tilde f`` over the sample (occupied bins, weighted by counts). The
# ``h_j`` drive the degree-0 adaptive KDE; ``g`` drives the degree-≥1 balloon bandwidth. O(nbins²).
function _abramson(counts, centers, n, h, kernel)
  nb = length(centers)
  pilot = zeros(Float64, nb)
  @inbounds for j in 1:nb
    iszero(counts[j]) && continue
    pilot[j] = _kde_h(counts, centers, n, centers[j], h, kernel)
  end
  logg = 0.0
  @inbounds for j in 1:nb
    (counts[j] > 0 && pilot[j] > 0) || continue
    logg += counts[j] * log(pilot[j])
  end
  g = exp(logg / n)
  hj = [pilot[j] > 0 ? h * sqrt(g / pilot[j]) : h for j in 1:nb]
  return hj, g
end


# Degree-0 adaptive KDE at x with per-bin (sample-point) bandwidths ``h_j``.
function _kde_adaptive(counts, centers, hj, n, x, kernel)
  acc = 0.0
  @inbounds for j in eachindex(centers)
    cj = counts[j]
    iszero(cj) && continue
    acc += cj * kernel((x - centers[j]) / hj[j]) / hj[j]
  end
  return acc / n
end


# Degree ≥ 1 — local-log-polynomial density by binned **local Poisson likelihood** (Loader 1996), at a
# single point-wise (balloon) bandwidth ``h_x``. Model the bin counts ``n_j`` as Poisson with
# ``\log \mathbb{E}[n_j] = \log(nb) + \sum_{l=0}^{d}\beta_l (c_j-x)^l`` and weights ``K((c_j-x)/h_x)``;
# Newton/IRLS (concave) gives ``\hat f(x) = e^{\beta_0}``. Warm-started from the KDE value `f0`; falls
# back to it where the local design can't support degree d (singular info matrix) or the fit diverges.
function _local_log(counts, centers, hx, n, b, x, degree, kernel, f0)
  (f0 > 0 && isfinite(f0)) || return f0
  invhx = inv(hx)
  p = degree + 1
  β = zeros(p)
  β[1] = log(f0)
  nb_off = n * b
  grad = zeros(p)
  H = zeros(p, p)
  pw = zeros(p)
  for _ in 1:25
    fill!(grad, 0.0)
    fill!(H, 0.0)
    @inbounds for j in eachindex(centers)
      d = centers[j] - x
      w = kernel(d * invhx)
      w < 1e-10 && continue
      pw[1] = 1.0
      η = β[1]
      acc = 1.0
      for l in 2:p
        acc *= d
        pw[l] = acc
        η += β[l] * acc
      end
      μ = nb_off * exp(η)
      r = counts[j] - μ
      for l in 1:p
        grad[l] += w * r * pw[l]
        wμpl = w * μ * pw[l]
        for m in 1:p
          H[l, m] += wμpl * pw[m]
        end
      end
    end
    δ = try
      H \ grad
    catch
      return f0
    end
    all(isfinite, δ) || return f0
    β .+= δ
    maximum(abs, δ) < 1e-8 && break
  end
  fx = exp(β[1])
  # a bounded refinement of the KDE: a fit that departs by more than ×5 is an extrapolation artifact
  # (a log-polynomial run wild in a sparse / wide-window tail), so fall back to the KDE there.
  return (isfinite(fx) && f0 / 5 ≤ fx ≤ 5 * f0) ? fx : f0
end


# The point-wise (balloon) bandwidth at x for the local-log fit: fixed `h`, or the Abramson factor
# ``h (\tilde f(x)/g)^{-1/2}`` evaluated *at the query point* (so the fit window stays local). The factor
# is clipped at 5 so the window does not blow up in the deep tail (where ``\tilde f(x)\to 0``) and the
# local polynomial extrapolate upward.
@inline _balloon(h, g, f0, adaptive) = (adaptive && f0 > 0) ? h * min(sqrt(g / f0), 5.0) : float(h)


"""
    density(summary::Histogram1D, x; degree = 0, h = nothing, adaptive = false, kernel = gaussian_kernel)   # ≡ f_X

**Estimand.** the marginal (Lebesgue) density ``f_X = \\mathrm{d}\\mu_X/\\mathrm{d}\\lambda`` of the
conditioning variable ``X``, defined by ``P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda``.

**Estimator.** a binned local smoother of the counts ``n_j`` at centers ``c_j``, with two orthogonal
knobs:

- **`degree`** — the local-polynomial order. `degree = 0` (default) is the binned KDE / WARPing
  ``\\hat f(x) = \\frac1n \\sum_j (n_j/h_j)\\, K((x-c_j)/h_j)`` (``\\int K = 1``). `degree ≥ 1` is the
  **local-log-polynomial** (local-likelihood) estimator (Loader 1996): fit ``\\log f`` by a
  degree-`degree` polynomial through a binned **local Poisson** likelihood and report
  ``\\hat f(x) = e^{\\hat\\beta_0}``. `degree = 1` (log-linear) is **exact for an exponential density**,
  so it tracks the exponential / stretched-exponential **tails** a flat KDE smears; `degree = 2`
  (log-quadratic) is exact for a Gaussian. It falls back to the KDE where the local fit is degenerate or
  runs away in a sparse tail.
- **`adaptive`** — the bandwidth. `false` (default) uses the fixed `h`; `true` uses the **Abramson**
  variable bandwidth — the kernel width scales by ``(\\tilde f/g)^{-1/2}`` (a pilot ``\\tilde f`` at `h`,
  ``g`` its geometric mean), widening in the sparse tails: far less tail variance (Abramson 1982). `h`
  is the fixed bandwidth or the pilot; it **defaults to the global Silverman rule**. (The KDE uses the
  sample-point ``h_j`` per bin; the local-log fit uses one balloon ``h(x)`` at the query.)

The four corners (fixed/variable ``\\times`` KDE/local-log) and everything between are one call. `kernel`
is any ``K`` with ``\\int K = 1`` ([`gaussian_kernel`](@ref GBGStreamingStatistics.DensityEstimation.gaussian_kernel)
default, [`epanechnikov_kernel`](@ref GBGStreamingStatistics.DensityEstimation.epanechnikov_kernel)); a
vector of query points returns a vector. The FFT amortization of the `degree = 0`, fixed-`h` convolution
over a bandwidth grid is forthcoming (`CLAUDE.md` roadmap).

```jldoctest
julia> h = Histogram1D(-1.0, 1.0, 2);

julia> observe!(h, [0.4]);          # one sample, binned at center 0.5

julia> round(density(h, 0.5; h = 0.2); digits = 4)
1.9947
```

References: I. S. Abramson (1982), [*On bandwidth variation in kernel estimates—a square root
law*](https://doi.org/10.1214/aos/1176345986), Ann. Statist. 10(4) (variable bandwidth);
C. R. Loader (1996), [*Local likelihood density estimation*](https://doi.org/10.1214/aos/1032298287),
Ann. Statist. 24(4) (the local-log-polynomial estimator).
"""
function density(summary::Histogram1D, x::Real; degree::Integer = 0, h = nothing,
                 adaptive::Bool = false, kernel = gaussian_kernel)
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  n = _total(summary)
  iszero(n) && return 0.0
  hh = h === nothing ? _silverman(summary) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = [_bincenter(summary, j) for j in 1:length(counts)]
  if degree == 0 && !adaptive
    return _kde_h(counts, centers, n, x, hh, kernel)
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    return _kde_adaptive(counts, centers, hj, n, x, kernel)
  end
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  f0 = _kde_h(counts, centers, n, x, hh, kernel)
  return _local_log(counts, centers, _balloon(hh, g, f0, adaptive), n, inv(summary.scale),
                    x, degree, kernel, f0)
end

function density(summary::Histogram1D, xs::AbstractVector; degree::Integer = 0, h = nothing,
                 adaptive::Bool = false, kernel = gaussian_kernel)
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  n = _total(summary)
  iszero(n) && return zeros(Float64, length(xs))
  hh = h === nothing ? _silverman(summary) : float(h)
  hh > 0 || throw(ArgumentError("bandwidth h must be positive (got h = $h)"))
  counts = Array(summary.counts)
  centers = [_bincenter(summary, j) for j in 1:length(counts)]
  if degree == 0 && !adaptive
    return [_kde_h(counts, centers, n, x, hh, kernel) for x in xs]
  elseif degree == 0
    hj, _ = _abramson(counts, centers, n, hh, kernel)
    return [_kde_adaptive(counts, centers, hj, n, x, kernel) for x in xs]
  end
  b = inv(summary.scale)
  g = adaptive ? _abramson(counts, centers, n, hh, kernel)[2] : 1.0
  return [begin
            f0 = _kde_h(counts, centers, n, x, hh, kernel)
            _local_log(counts, centers, _balloon(hh, g, f0, adaptive), n, b, x, degree, kernel, f0)
          end for x in xs]
end

end # module DensityEstimation
