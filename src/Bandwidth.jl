# Bandwidth — §4.5: bandwidth selection (the kernel/degree normal-reference, Sheather–Jones — classic
# and improved (DCT) — fit-based AIC, the Chen asymmetric-kernel plug-in, and the response-aware local
# plug-in (rule = :local) now; the FFT-amortized smoothing is forthcoming). The smoothing is a convolution of the binned summary with the kernel; the
# FFT computes it in O(B log B) and amortizes over bandwidth (`docs/ARCHITECTURE.md` §4.5) — hence the
# FFTW dependency, used once that lands.

module Bandwidth

using AbstractFFTs
using FFTW

using ..Histograms: Histogram1D, Histogram2D, _total, _sample_std, _bincenter, _column_sums
using ..LocalPolynomial: local_fit
using ..Kernels: gaussian_kernel, _quad
using ..Regression: _residual_response
using ..MarginalDensity: _normal_reference_h, _kde_h, AsymmetricKernel, _asymmetric_default_b

export bandwidth


# --- Sheather–Jones solve-the-equation plug-in (`rule = :sj`), from the binned summary ---
# The data-driven second-order (degree 0/1) bandwidth: estimate the density curvature ψ₄ = R(f'') and
# its ψ₆ pilot from the bin counts as Hermite-weighted reductions, then solve the Sheather–Jones
# equation (Sheather & Jones 1991), mirroring base-R `bw.SJ(method = "ste")`. The functional pilots are
# Gaussian (as in base-R / locfit); the kernel enters only through c₁ = R(K)/(μ₂(K)²n), so the result
# stays kernel-aware. Both functionals are reductions of the bin counts, hence mergeable.

# Gaussian-derivative weights φ⁽⁴⁾, φ⁽⁶⁾ (Hermite × Gaussian), for the ψ₄, ψ₆ functionals.
_phi4(u) = (u^4 - 6u^2 + 3) * exp(-u^2 / 2) / sqrt(2π)
_phi6(u) = (u^6 - 15u^4 + 45u^2 - 15) * exp(-u^2 / 2) / sqrt(2π)


# Distance histogram of the bin counts: dist[m+1] = Σ_j counts_j·counts_{j+m} (the autocorrelation — the
# pair count at bin-lag m). Computed once; each functional evaluation is then O(nbins). O(nbins²) here;
# the FFT autocorrelation is the §4.5 performance lever, forthcoming.
function _distance_hist(counts)
  nbins = length(counts)
  dist = zeros(Float64, nbins)
  @inbounds for j in 1:nbins
    cj = counts[j]
    cj == 0 && continue
    for k in j:nbins
      dist[k - j + 1] += cj * counts[k]
    end
  end
  return dist
end


# The density functional ψ_r(g) = (1/(n²g^{r+1}))[Σ_{j,k} counts_j counts_k φ⁽ʳ⁾((c_j−c_k)/g) − nφ⁽ʳ⁾(0)],
# the n self-pairs removed, from the distance histogram. φ⁽ʳ⁾ is even, so bin-lags ±m fold to 2·dist[m].
function _psi(dist, binwidth, n, g, phir, r)
  acc = dist[1] * phir(0.0)
  @inbounds for m in 1:length(dist)-1
    d = dist[m + 1]
    d == 0 && continue
    acc += 2 * d * phir(m * binwidth / g)
  end
  acc -= n * phir(0.0)
  return acc / (n^2 * g^(r + 1))
end


# Binned empirical p-quantile (linear interpolation in the crossing bin) — for the robust SJ scale.
function _binned_quantile(summary::Histogram1D, p)
  counts = summary.counts
  target = p * _total(summary)
  cum = zero(eltype(counts))
  @inbounds for j in eachindex(counts)
    cum += counts[j]
    if cum ≥ target
      frac = counts[j] > 0 ? (target - (cum - counts[j])) / counts[j] : 0.5
      return summary.lower + (j - 1 + frac) / summary.scale
    end
  end
  return summary.upper
end


# canonical bandwidth rescale δ(K)/δ(φ) = [2√π·R(K)/μ₂(K)²]^(1/5) (= 1 for the Gaussian): the factor by
# which any eval kernel's AMISE-optimal bandwidth differs from the Gaussian's. Shared by the data-driven
# selectors, which estimate the data's curvature with a Gaussian and apply the kernel as this rescale.
_canonical_ratio(kernel) =
  (2 * sqrt(π) * _quad(v -> kernel(v)^2) / _quad(v -> v^2 * kernel(v))^2)^(1 / 5)


function _sheather_jones(summary::Histogram1D, kernel)
  counts = Float64.(Array(summary.counts))
  n = float(_total(summary))
  binwidth = inv(summary.scale)
  σ = _sample_std(summary)
  iqr = _binned_quantile(summary, 0.75) - _binned_quantile(summary, 0.25)
  scale = (iqr > 0) ? min(σ, iqr / 1.349) : σ
  dist = _distance_hist(counts)
  c1 = inv(2 * sqrt(π) * n)             # standard (Gaussian-pilot) SJ: R(φ)/(μ₂(φ)²n) = 1/(2√πn)
  a = 1.24 * scale * n^(-1 / 7)
  b = 1.23 * scale * n^(-1 / 9)
  TD = -_psi(dist, binwidth, n, b, _phi6, 6)
  SDa = _psi(dist, binwidth, n, a, _phi4, 4)
  (TD > 0 && SDa > 0) ||
    throw(ArgumentError("sample too sparse for Sheather–Jones; use rule = :normal"))
  alph2 = 1.357 * (SDa / TD)^(1 / 7)
  function fSD(h)                        # the Sheather–Jones equation residual
    ψ4 = _psi(dist, binwidth, n, alph2 * h^(5 / 7), _phi4, 4)
    ψ4 > 0 ? (c1 / ψ4)^(1 / 5) - h : Inf   # ψ₄ ≤ 0: pilot too small ⇒ push the root toward larger h
  end
  # bracket around the normal-reference scale, expand if needed, then bisect to the root
  lo = 0.1 * scale * n^(-1 / 5)
  hi = 2.0 * scale * n^(-1 / 5)
  flo, fhi = fSD(lo), fSD(hi)
  itry = 0
  while flo * fhi > 0 && itry < 50
    if isodd(itry)
      hi *= 1.2; fhi = fSD(hi)
    else
      lo /= 1.2; flo = fSD(lo)
    end
    itry += 1
  end
  flo * fhi ≤ 0 ||
    throw(ArgumentError("no Sheather–Jones solution bracketed; use rule = :normal"))
  for _ in 1:60
    mid = 0.5 * (lo + hi)
    fm = fSD(mid)
    if flo * fm ≤ 0
      hi = mid; fhi = fm
    else
      lo = mid; flo = fm
    end
  end
  hG = 0.5 * (lo + hi)
  return hG * _canonical_ratio(kernel)   # the eval kernel enters as a clean canonical δ(K) rescale
end


# --- fit-based selector over an h-grid (rule = :aic), from the binned summary (degree 0/1) ---
# Minimise the AIC of the Gaussian KDE over a log-spaced h-grid around the normal reference, read off the
# bin counts; the eval kernel then enters as the canonical δ(K) rescale. The criterion is a deterministic
# functional of the counts, hence mergeable. The cross-validation selectors (likelihood LCV and
# least-squares LSCV) are **deliberately absent**: both are leave-one-out, and binning ties (many samples
# share a bin) make the same-bin term drive the CV optimum to h → 0 — a recorded coverage gap (a within-
# bin-spread correction from the power-sums could rescue them later). AIC needs no leave-one-out — both
# its log-likelihood and df-penalty diverge as h → 0 — so it has a stable interior optimum; `:sj` is the
# other reliable data-driven selector.

_gauss(u) = exp(-u^2 / 2) / sqrt(2π)


# AIC of the Gaussian KDE: −2 Σ_j n_j log f̂(c_j) + 2·df, df = (φ(0)/h) Σ_j n_j/(n f̂(c_j)) (the KDE
# self-influence trace). Both terms diverge as h → 0, so the minimiser is a stable interior optimum. The
# KDE convolution is summed only over bins within the kernel reach (±5h), so it is O(occupied · h⁻¹).
function _aic_gauss(counts, occ, n, h, binwidth)
  invh = inv(h)
  reach = ceil(Int, 5 * h / binwidth)
  nb = length(counts)
  ll = 0.0
  df = 0.0
  @inbounds for j in occ
    fj = 0.0
    for k in max(1, j - reach):min(nb, j + reach)
      ck = counts[k]
      ck == 0 && continue
      fj += ck * _gauss((j - k) * binwidth * invh)
    end
    fj *= invh / n
    fj > 0 || continue
    ll += counts[j] * log(fj)
    df += counts[j] / (n * fj)
  end
  return -2 * ll + 2 * _gauss(0.0) * invh * df
end


function _aic_select(summary::Histogram1D, kernel)
  counts = Float64.(Array(summary.counts))
  n = float(_total(summary))
  binwidth = inv(summary.scale)
  bw0 = _normal_reference_h(summary, gaussian_kernel, 0)   # the Gaussian normal reference centres the grid
  occ = findall(>(0), counts)
  best_h = bw0
  best_val = Inf
  for h in bw0 .* exp.(range(log(0.2), log(3.0); length = 50))
    val = _aic_gauss(counts, occ, n, h, binwidth)
    if val < best_val
      best_val = val
      best_h = h
    end
  end
  return best_h * _canonical_ratio(kernel)
end


# --- Improved Sheather–Jones (Botev 2010), `rule = :isj`, from the binned summary (degree 0/1) ---
# The solve-the-equation bandwidth read off the DCT of the bin counts: the curvature norms ‖f^(s)‖² are
# evaluated from the cosine coefficients (Plancherel) and the AMISE fixed point t = ξ·γ⁵(t) is solved by
# root-finding; h = √t·(upper−lower). Unlike the normal-reference pilot of `:sj`, the curvature is
# estimated nonparametrically, so ISJ stays robust on strongly multimodal data (Botev, Grotowski & Kroese 2010). A
# deterministic functional of the counts (mergeable); the eval kernel enters as the δ(K) rescale. The
# I²^ℓ·exp(−I²π²t) terms are formed in log-space, since I²^7 overflows Float64 at fine grids.
function _isj_fixed_point(t, N, Isq, logIsq, a2)
  ell = 7                                                # 5 plug-in stages (Botev's recommendation)
  f = 0.5 * π^(2ell) * sum(@. exp(ell * logIsq - Isq * π^2 * t) * a2)
  f ≤ 0 && return -one(f)
  for s in (ell-1):-1:2
    K0 = prod(1:2:(2s-1)) / sqrt(2π)                     # (2s−1)!! / √(2π) — the Gaussian-derivative norm
    c = (1 + (0.5)^(s + 0.5)) / 3
    time = (2 * c * K0 / (N * f))^(2 / (3 + 2s))
    f = 0.5 * π^(2s) * sum(@. exp(s * logIsq - Isq * π^2 * time) * a2)
  end
  return t - (2 * N * sqrt(π) * f)^(-2 / 5)
end


function _isj(summary::Histogram1D, kernel)
  counts = Float64.(Array(summary.counts))
  nbins = length(counts)
  N = float(_total(summary))
  a = FFTW.r2r(counts ./ N, FFTW.REDFT10)               # type-2 DCT of the binned density (sums to 1)
  a2 = a[2:end] .^ 2                                     # squared coefficients, DC term dropped
  Isq = Float64[k^2 for k in 1:nbins-1]
  logIsq = log.(Isq)
  fp(t) = _isj_fixed_point(t, N, Isq, logIsq, a2)
  tol = 1e-11 + 0.01 * (clamp(N, 50.0, 1050.0) - 50) / 1000
  lo, hi = 0.0, tol
  flo, fhi = fp(lo), fp(hi)
  itry = 0
  while flo * fhi > 0 && itry < 60                       # expand the bracket until it straddles the root
    hi *= 2; fhi = fp(hi); itry += 1
  end
  flo * fhi ≤ 0 ||
    throw(ArgumentError("Improved Sheather–Jones did not bracket a root (too sparse); use rule = :sj"))
  for _ in 1:100
    mid = 0.5 * (lo + hi)
    fm = fp(mid)
    flo * fm ≤ 0 ? (hi = mid; fhi = fm) : (lo = mid; flo = fm)
  end
  return sqrt(0.5 * (lo + hi)) * (summary.upper - summary.lower) * _canonical_ratio(kernel)
end


# --- Chen plug-in for the asymmetric (gamma / beta) kernels (rule = :plugin), from the binned summary ---
# The asymmetric kernels take a single scale `b`, not the symmetric `h`. Their AMISE (Chen 1999, 2000;
# modified f̂₂ shapes) is A·b² + C·b^(−1/2)/n with A = ¼∫w(x)²f''(x)² and C = (1/(2√π))𝔼[w(X)^(−1/2)],
# where w is the kernel's AMISE weight (x for gamma, x(1−x) for beta — the leading bias is ½b·w·f'' and
# the variance ∝ w^(−1/2)). Minimising gives the analog of Silverman's plug-in,
#   b⋆ = [C/(4An)]^(2/5) = [𝔼[w(X)^(−1/2)] / (2√π·n·∫w²f'')]^(2/5).
# 𝔼[w(X)^(−1/2)] is an exact moment of the counts; the curvature ∫w²f'' is the data-driven part.

# ∫ w(x)² f''(x)² dx by grid quadrature of the pilot KDE's second derivative f̂''_g (Gaussian pilot,
# φ''(u) = (u²−1)φ(u)), with the variance inflation R(φ'')·𝔼[w²]/(n g⁵) removed (R(φ'') = 3/(8√π)). On a
# bounded / positive support this is the robust route: w(x) → 0 at the boundary suppresses the symmetric
# pilot's edge leakage, and — unlike the φ⁴ integration-by-parts U-statistic of `_psi` — it assumes **no**
# boundary-vanishing, so it stays valid where f(0) > 0 (e.g. the exponential, where the U-statistic fails).
function _curvature_grid(counts, centers, occ, n, g, weight, lower, upper)
  npts = 2048
  step = (upper - lower) / (npts - 1)
  inv_g = inv(g)
  pref = inv(n * g^3) / sqrt(2π)
  acc = 0.0
  @inbounds for m in 0:npts-1
    x = lower + m * step
    fpp = 0.0
    for j in occ
      u = (x - centers[j]) * inv_g
      fpp += counts[j] * (u * u - 1) * exp(-u * u / 2)
    end
    fpp *= pref
    w = weight(x)
    acc += w * w * fpp * fpp * step
  end
  ew2 = 0.0
  @inbounds for j in occ
    ew2 += counts[j] * weight(centers[j])^2
  end
  bias = (3 / (8 * sqrt(π))) * (ew2 / n) / (n * g^5)   # the f̂'' variance inflation, R(φ'')·𝔼[w²]/(n g⁵)
  return max(acc - bias, 0.0)
end


function _plugin_asymmetric(summary::Histogram1D, kernel::AsymmetricKernel)
  counts = Float64.(Array(summary.counts))
  n = float(_total(summary))
  centers = [_bincenter(summary, j) for j in eachindex(counts)]
  klower, kupper = kernel.support
  occ = [j for j in eachindex(counts) if counts[j] > 0 && klower < centers[j] < kupper]
  isempty(occ) &&
    throw(ArgumentError("no occupied in-support bins for the plug-in bandwidth; use rule = :normal"))
  weight = kernel.weight
  weight === nothing &&
    throw(ArgumentError("this `AsymmetricKernel` was built without an AMISE weight, so `rule = :plugin` has no `w(x)` to use; supply it as `AsymmetricKernel(pdf, support, weight)`, or select `rule = :normal`"))
  σ = _sample_std(summary)
  iqr = _binned_quantile(summary, 0.75) - _binned_quantile(summary, 0.25)
  scale = (iqr > 0) ? min(σ, iqr / 1.349) : σ
  # the f''-estimation pilot: a normal-reference scale at the second-derivative rate n^(−1/9). The
  # constant (≈ 1.3, just above the Gaussian-reference value (4/7)^{1/9} ≈ 0.94) is set — as Silverman's
  # 1.06 is — so the bias-corrected curvature gives a bandwidth that is robust across positive densities
  # (never materially worse than the rule of thumb, ≈ 1.4–1.7× lower ISE on skewed data). The
  # fully-adaptive solve-the-equation pilot is a recorded refinement (CLAUDE.md roadmap).
  g = 1.3 * scale * n^(-1 / 9)
  lower = max(klower, float(summary.lower))
  upper = min(kupper, float(summary.upper))
  curvature = _curvature_grid(counts, centers, occ, n, g, weight, lower, upper)
  curvature > 0 ||
    throw(ArgumentError("plug-in curvature estimate ≤ 0 (sample too sparse or nearly flat); use rule = :normal"))
  ew = sum(counts[j] * weight(centers[j])^(-1 / 2) for j in occ) / n   # 𝔼[w(X)^(−1/2)]
  return (ew / (2 * sqrt(π) * n * curvature))^(2 / 5)
end


"""
    bandwidth(summary::Histogram1D; rule = :normal, kernel = gaussian_kernel,
              degree = 0) -> h

**Estimand.** the AMISE-optimal smoothing bandwidth for the marginal density
[`density`](@ref GBGStreamingStatistics.MarginalDensity.density) of `summary` — the ``h`` minimising
``\\mathrm{MISE}(\\hat f_h) = \\int \\mathbb{E}\\bigl[(\\hat f_h(x) - f_X(x))^2\\bigr]\\,\\mathrm{d}x``
for the smoothing `kernel` ``K`` and local-polynomial `degree` ``d``.

**Estimator** (`rule = :normal`, the kernel- and degree-aware normal-reference plug-in). From the
recovered ``\\hat\\sigma`` and the sample size ``n``,

```math
h = C(K, d)\\,\\hat\\sigma\\,n^{-1/(2p+1)},
```

where ``p`` is the estimator's leading **bias order** — ``p = 2`` for `degree` 0 or 1 (the rate
``n^{-1/5}``), ``p = 4`` for `degree` 2 or 3 (a higher-order, ``O(h^4)``, fit, hence the wider, more
slowly shrinking rate ``n^{-1/9}``). At `degree` 0/1 the constant is the **canonical bandwidth**

```math
\\begin{aligned}
C(K, d) &= (8\\sqrt{\\pi}/3)^{1/5}\\,\\delta(K), \\\\
\\delta(K) &= \\bigl[R(K)/\\mu_2(K)^2\\bigr]^{1/5}, \\\\
R(K) &= \\int K^2, \\qquad \\mu_2(K) = \\int u^2 K(u)\\,\\mathrm{d}u,
\\end{aligned}
```

and at `degree` 2/3 the analogous fourth-order **equivalent-kernel** constant (Loader 1996). ``R(K)``
and ``\\mu_2(K)`` are computed by quadrature, so any ``K`` is supported. For the Gaussian kernel at
`degree` 0 this is **Silverman's rule** ``(4/3)^{1/5}\\hat\\sigma\\,n^{-1/5}``; for the Epanechnikov
kernel it is ``\\approx 2.34\\,\\hat\\sigma\\,n^{-1/5}`` (``2.21\\times`` the Gaussian ``h``). Every
input is read back from the *one* stored summary, post hoc.

**Estimator** (`rule = :sj`, the data-driven Sheather–Jones plug-in; `degree` 0/1). The curvature
``\\psi_4 = R(f'')`` and its ``\\psi_6`` pilot are estimated from the bin counts as Hermite-weighted
reductions, and the **solve-the-equation** bandwidth (Sheather & Jones 1991) is found by root-finding;
the eval kernel then enters as the same canonical ``\\delta(K)`` rescale as for `:normal` (the curvature
is a Gaussian-pilot estimate, kernel-independent). It adapts to the true shape — on a bimodal sample,
where the normal reference oversmooths from an inflated ``\\hat\\sigma``, it picks a much smaller ``h``.
Both functionals are reductions of the bin counts, so they are mergeable; a `degree` ``\\ge 2`` request
is a recorded coverage gap (use `rule = :normal`).

**Estimator** (`rule = :isj`, the **Improved Sheather–Jones** plug-in; Botev, Grotowski & Kroese 2010;
`degree` 0/1). The same solve-the-equation idea, but the curvature norms ``\\|f^{(s)}\\|^2`` are read off
the **DCT** of the bin counts (Plancherel) rather than from a normal-reference pilot — so it stays
reliable on strongly multimodal data, where the pilot scale of `:sj` can mislead. Mergeable (a functional
of the counts' cosine transform); the eval kernel again enters as the ``\\delta(K)`` rescale.

**Estimator** (`rule = :aic`, the fit-based information criterion; `degree` 0/1). Over a log-spaced
``h``-grid around the normal reference, minimise ``\\mathrm{AIC}(h) = -2\\sum_j n_j \\log \\hat f_h(c_j) + 2\\,\\mathrm{df}(h)`` of the binned **Gaussian** KDE, with ``\\mathrm{df}(h) = (\\varphi(0)/h)\\sum_j n_j/(n\\hat f_h(c_j))`` the kernel self-influence trace (``\\varphi`` the Gaussian). Both terms diverge as ``h \\to 0``, so the optimum is a stable interior point; like `:sj` it is a deterministic functional of the counts (mergeable), and the eval kernel enters only as the canonical ``\\delta(K)`` rescale of the selected ``h``. The leave-one-out **cross-validation** selectors (likelihood-CV and least-squares LSCV) are a deliberate coverage gap — binning ties drive their optimum to ``h \\to 0`` (use `:sj` for a data-driven alternative).

**Asymmetric kernels** ([`gamma_kernel`](@ref GBGStreamingStatistics.MarginalDensity.gamma_kernel),
[`beta_kernel`](@ref GBGStreamingStatistics.MarginalDensity.beta_kernel)) carry their width in the
varying shape, so they take their own scale ``b`` (rate ``n^{-2/5}``), not the symmetric ``h``, and accept
two rules. `rule = :normal` (default) is the ``\\hat\\sigma\\,n^{-2/5}`` rule of thumb. `rule = :plugin`
is the **Chen plug-in**, the asymmetric analog of Silverman's rule: from the AMISE ``A b^2 + C b^{-1/2}/n``
(Chen 1999, 2000; the modified ``f_2`` kernels) with ``A = \\tfrac14\\int w(x)^2 f''(x)^2\\,\\mathrm{d}x``
and ``C = (2\\sqrt\\pi)^{-1}\\,\\mathbb{E}[w(X)^{-1/2}]``, where ``w(x) = x`` (gamma) or ``x(1-x)`` (beta)
is the kernel's local-scale weight,

```math
b^\\star = \\left[\\frac{\\mathbb{E}[w(X)^{-1/2}]}{2\\sqrt\\pi\\,n \\int w(x)^2 f''(x)^2\\,\\mathrm{d}x}\\right]^{2/5}.
```

``\\mathbb{E}[w(X)^{-1/2}]`` is an exact moment of the bin counts; the curvature ``\\int w^2 (f'')^2`` is
estimated from a pilot KDE's second derivative by grid quadrature (valid at the boundary because
``w \\to 0`` there, unlike an integration-by-parts functional). It lowers the ISE by ``\\approx
1.4\\text{–}1.7\\times`` on skewed positive data — where the rule of thumb mis-scales — and is otherwise near
it. The symmetric selectors `:sj` / `:isj` / `:aic` do not apply to an asymmetric kernel; a fully
data-adaptive curvature pilot is a recorded refinement.

The response-aware *local* plug-in ``h^\\star(x) = \\bigl[R(K)\\,\\sigma^2(x) /
(\\mu_2(K)^2\\, m''(x)^2\\, f_X(x)\\, n)\\bigr]^{1/5}`` — each conditioning value its own bandwidth from
the conditional variance ``\\sigma^2(x)`` and regression curvature ``m''(x)`` — is the separate
joint-summary selector `bandwidth(summary::Histogram2D, x; rule = :local)` (below).

References:

- B. W. Silverman (1986), [*Density Estimation for Statistics and Data Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall (the normal-reference rule and the canonical bandwidth).
- S. P. Sheather & M. C. Jones (1991), [*A reliable data-based bandwidth selection method for kernel density estimation*](https://doi.org/10.1111/j.2517-6161.1991.tb01857.x), J. R. Statist. Soc. B 53(3) (the solve-the-equation plug-in).
- Z. I. Botev, J. F. Grotowski & D. P. Kroese (2010), [*Kernel density estimation via diffusion*](https://doi.org/10.1214/10-AOS799), Ann. Statist. 38(5) (the DCT solve-the-equation, `:isj`).
- C. R. Loader (1996), [*Local likelihood density estimation*](https://doi.org/10.1214/aos/1032298287), Ann. Statist. 24(4) (the degree / equivalent-kernel factor).
"""
function bandwidth(summary::Histogram1D; rule::Symbol = :normal,
                   kernel = gaussian_kernel, degree::Integer = 0)
  n = _total(summary)
  n ≥ 2 || throw(ArgumentError("bandwidth needs at least 2 samples (got n = $(Int(n)))"))
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  if kernel isa AsymmetricKernel
    degree == 0 ||
      throw(ArgumentError("the asymmetric kernel is a degree-0 estimator; `bandwidth(…; kernel = gamma_kernel)` takes no `degree` (got degree = $degree)"))
    rule === :normal && return _asymmetric_default_b(summary)
    rule === :plugin && return _plugin_asymmetric(summary, kernel)
    throw(ArgumentError("for an asymmetric kernel use rule = :normal (the σ̂ n^(−2/5) rule of thumb) " *
                        "or :plugin (the Chen curvature plug-in); the symmetric-kernel selectors " *
                        ":sj / :isj / :aic do not apply (got :$rule)"))
  end
  if rule === :normal
    return _normal_reference_h(summary, kernel, degree)
  elseif rule === :sj
    (degree == 0 || degree == 1) ||
      throw(ArgumentError("Sheather–Jones (rule = :sj) is a second-order selector (degree ≤ 1); " *
                          "for degree ≥ 2 use rule = :normal (got degree = $degree)"))
    return _sheather_jones(summary, kernel)
  elseif rule === :isj
    (degree == 0 || degree == 1) ||
      throw(ArgumentError("Improved Sheather–Jones (rule = :isj) is a degree-0/1 selector; " *
                          "for degree ≥ 2 use rule = :normal (got degree = $degree)"))
    return _isj(summary, kernel)
  elseif rule === :aic
    (degree == 0 || degree == 1) ||
      throw(ArgumentError("the fit-based selector rule = :aic is a degree-0/1 (KDE) selector; " *
                          "for degree ≥ 2 use rule = :normal (got degree = $degree)"))
    return _aic_select(summary, kernel)
  end
  throw(ArgumentError("unknown bandwidth rule :$rule (expected :normal, :sj, :isj, or :aic)"))
end


# --- response-aware local plug-in (rule = :local), from the joint summary (Histogram2D) ---
# Balancing the local-linear squared bias (∝ h⁴ m''(x)²) against the variance (∝ σ²(x)/(n h f_X(x)))
# gives the AMISE-optimal *local* bandwidth h⋆(x) = [R(K) σ²(x) / (μ₂(K)² m''(x)² f_X(x) n)]^(1/5). The
# three x-local pilots — σ²(x), m''(x) = 2β₂ (a degree-2 fit), f_X(x) (the X-marginal KDE) — are read
# from the one joint summary at a single pilot bandwidth, wider than h⋆ (the n^(−1/9) curvature rate, so
# the second derivative is stable). Where the regression is locally linear (m'' → 0) the optimal h is
# unbounded, so it is capped at the X-range.

# The default pilot: σ̂_X n^(−1/9), the rule-of-thumb scale at the second-derivative (curvature) rate —
# wider than the final h⋆ ∝ n^(−1/5), since estimating m'' needs a more stable (wider) window. σ̂_X is
# the binned X-marginal standard deviation (the column centers, adequate for a pilot).
function _local_pilot(cx, N, n)
  xbar = 0.0
  @inbounds for j in eachindex(cx)
    xbar += N[j] * cx[j]
  end
  xbar /= n
  v = 0.0
  @inbounds for j in eachindex(cx)
    v += N[j] * (cx[j] - xbar)^2
  end
  return sqrt(v / n) * n^(-1 / 9)
end


function _hstar_at(cx, N, sumY, SR2, n, x, kernel, hpilot, RK, mu2, xrange)
  fX = _kde_h(N, cx, n, x, hpilot, kernel)
  fX > 0 || return oftype(fX, NaN)                       # x outside the data's support
  sigma2 = local_fit(cx, N, SR2, x, hpilot, kernel, 1)[1]      # Fan–Yao conditional variance σ²(x)
  mpp = 2 * local_fit(cx, N, sumY, x, hpilot, kernel, 2)[3]    # m''(x) = 2 β₂
  (isfinite(sigma2) && isfinite(mpp)) || return oftype(fX, NaN)
  sigma2 = max(sigma2, 0.0)
  denom = mu2^2 * mpp^2 * fX * n
  denom > 0 || return oftype(fX, xrange)                 # locally linear (m'' ≈ 0) ⇒ h unbounded; cap
  return min((RK * sigma2 / denom)^(1 / 5), xrange)
end


function _local_setup(summary::Histogram2D, kernel, hpilot)
  cx, N, sumY, sumY2 = _column_sums(summary)
  n = sum(N)
  n ≥ 2 || throw(ArgumentError("bandwidth needs at least 2 samples (got n = $n)"))
  hp = hpilot === nothing ? _local_pilot(cx, N, n) : float(hpilot)
  hp > 0 || throw(ArgumentError("pilot bandwidth hpilot must be positive (got hpilot = $hp)"))
  RK = _quad(v -> kernel(v)^2)
  mu2 = _quad(v -> v^2 * kernel(v))
  xrange = summary.xupper - summary.xlower
  SR2 = _residual_response(cx, N, sumY, sumY2, hp, kernel, 1)   # Fan–Yao residual response for σ²(x)
  return cx, N, sumY, SR2, n, hp, RK, mu2, xrange
end


"""
    bandwidth(summary::Histogram2D, x;
              rule = :local, kernel = gaussian_kernel, hpilot = nothing)   # ≡ h⋆(x)

**Estimand.** the response-aware **local** AMISE-optimal bandwidth ``h^\\star(x)`` for the conditional
mean [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean) ``m(x)`` — the
bandwidth that, at each conditioning value ``x``, minimises the local mean squared error of the
local-linear estimate ``\\hat m_h(x)`` by balancing its squared bias against its variance:

```math
h^\\star(x) = \\left[\\frac{R(K)\\,\\sigma^2(x)}{\\mu_2(K)^2\\,m''(x)^2\\,f_X(x)\\,n}\\right]^{1/5},
```

with ``\\sigma^2(x)`` the conditional variance, ``m''(x)`` the regression curvature, ``f_X(x)`` the
marginal density of the conditioning variable ``X``, ``n`` the sample size, and the kernel constants
``R(K) = \\int K^2``, ``\\mu_2(K) = \\int u^2 K(u)\\,\\mathrm{d}u``.

**Estimator** (`rule = :local`). The three ``x``-local quantities are read from the one joint `summary`,
each at a single pilot bandwidth: the conditional variance ``\\sigma^2(x)`` from the bias-reduced Fan &
Yao `conditional_variance` (its default residual estimator), the curvature ``m''(x) = 2\\beta_2`` from the
degree-2 local-polynomial fit, and ``f_X(x)`` from the binned KDE of the X-marginal. ``R(K)`` and
``\\mu_2(K)`` are computed by quadrature, so any `kernel` works. Because estimating a second derivative
needs a wider window than the final fit, the default pilot is the curvature-rate rule of thumb
``\\hat\\sigma_X\\, n^{-1/9}``; where the regression is locally linear (``m''(x) \\to 0``) the optimal
bandwidth is unbounded, so it is capped at the X-range. A vector `x` returns a vector; `NaN` marks an
``x`` off the data's support.

Arguments:

- `summary` — the joint [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D) of
  ``(X, Y)``.
- `x` — the conditioning value(s) at which ``h^\\star`` is evaluated (a scalar, or a vector for the
  whole bandwidth profile).
- `rule` — the selector; only `:local` (the response-aware local plug-in) is defined for a joint
  summary.
- `kernel` — the smoothing kernel ``K`` (any callable with ``\\int K = 1``; default `gaussian_kernel`).
- `hpilot` — the pilot bandwidth for the ``\\sigma^2(x)``, ``m''(x)``, ``f_X(x)`` estimates; `nothing`
  (default) uses the curvature-rate rule of thumb ``\\hat\\sigma_X\\, n^{-1/9}``.

Reference: J. Fan & I. Gijbels (1996), [*Local Polynomial Modelling and Its Applications*](https://doi.org/10.1201/9780203748725),
Chapman & Hall.
"""
function bandwidth(summary::Histogram2D, x::Real; rule::Symbol = :local,
                   kernel = gaussian_kernel, hpilot = nothing)
  rule === :local ||
    throw(ArgumentError("a joint summary supports only rule = :local (the response-aware local " *
                        "plug-in); got :$rule"))
  cx, N, sumY, SR2, n, hp, RK, mu2, xrange = _local_setup(summary, kernel, hpilot)
  return _hstar_at(cx, N, sumY, SR2, n, x, kernel, hp, RK, mu2, xrange)
end

function bandwidth(summary::Histogram2D, xs::AbstractVector; rule::Symbol = :local,
                   kernel = gaussian_kernel, hpilot = nothing)
  rule === :local ||
    throw(ArgumentError("a joint summary supports only rule = :local (the response-aware local " *
                        "plug-in); got :$rule"))
  cx, N, sumY, SR2, n, hp, RK, mu2, xrange = _local_setup(summary, kernel, hpilot)
  return [_hstar_at(cx, N, sumY, SR2, n, x, kernel, hp, RK, mu2, xrange) for x in xs]
end

end # module Bandwidth
