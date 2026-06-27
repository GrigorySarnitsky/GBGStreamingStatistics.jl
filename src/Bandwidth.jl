# Bandwidth — §4.5: bandwidth selection (the kernel/degree normal-reference and Sheather–Jones rules
# now; the fit-based LCV/AIC selectors, the response-aware local plug-in, and the FFT-amortized
# smoothing are forthcoming). The smoothing is a convolution of the binned summary with the kernel; the
# FFT computes it in O(B log B) and amortizes over bandwidth (`docs/ARCHITECTURE.md` §4.5) — hence the
# FFTW dependency, used once that lands.

module Bandwidth

using AbstractFFTs
using FFTW

using ..Histograms: Histogram1D, _total, _sample_std
using ..DensityEstimation: gaussian_kernel, _normal_reference_h, _quad

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
  # the eval kernel enters as a clean canonical rescale δ(K)/δ(φ) = [2√π·R(K)/μ₂(K)²]^(1/5) (= 1 for the
  # Gaussian) — the curvature ψ₄ is a property of the data, estimated once by the Gaussian-pilot ste.
  return hG * (2 * sqrt(π) * _quad(v -> kernel(v)^2) / _quad(v -> v^2 * kernel(v))^2)^(1 / 5)
end


"""
    bandwidth(summary::Histogram1D; rule = :normal, kernel = gaussian_kernel, degree = 0) -> h

**Estimand.** the AMISE-optimal smoothing bandwidth for the marginal density
[`density`](@ref GBGStreamingStatistics.DensityEstimation.density) of `summary` — the ``h`` minimising
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
C(K, d) = (8\\sqrt{\\pi}/3)^{1/5}\\,\\delta(K),
\\qquad
\\delta(K) = \\bigl[R(K)/\\mu_2(K)^2\\bigr]^{1/5},
\\quad R(K) = \\int K^2,\\ \\mu_2(K) = \\int u^2 K(u)\\,\\mathrm{d}u,
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

The fit-based **likelihood-CV / AIC** selectors (`rule = :lcv`, `:aic`) are forthcoming (`CLAUDE.md`
roadmap), as is the response-aware *local* plug-in ``h^\\star(x) = \\bigl[R(K)\\,\\sigma^2(x) /
(\\mu_2(K)^2\\, m''(x)^2\\, f_X(x)\\, n)\\bigr]^{1/5}`` — each response channel its own bandwidth from
the conditional variance ``\\sigma^2(x)`` and regression curvature ``m''(x)``.

References: B. W. Silverman (1986), [*Density Estimation for Statistics and Data
Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall (the normal-reference rule and the
canonical bandwidth); S. P. Sheather & M. C. Jones (1991), [*A reliable data-based bandwidth selection
method for kernel density estimation*](https://doi.org/10.1111/j.2517-6161.1991.tb01857.x), J. R.
Statist. Soc. B 53(3) (the solve-the-equation plug-in); C. R. Loader (1996), [*Local likelihood density
estimation*](https://doi.org/10.1214/aos/1032298287), Ann. Statist. 24(4) (the degree / equivalent-
kernel factor).
"""
function bandwidth(summary::Histogram1D; rule::Symbol = :normal,
                   kernel = gaussian_kernel, degree::Integer = 0)
  n = _total(summary)
  n ≥ 2 || throw(ArgumentError("bandwidth needs at least 2 samples (got n = $(Int(n)))"))
  degree ≥ 0 || throw(ArgumentError("degree must be ≥ 0 (got $degree)"))
  if rule === :normal
    return _normal_reference_h(summary, kernel, degree)
  elseif rule === :sj
    (degree == 0 || degree == 1) ||
      throw(ArgumentError("Sheather–Jones (rule = :sj) is a second-order selector (degree ≤ 1); " *
                          "for degree ≥ 2 use rule = :normal (got degree = $degree)"))
    return _sheather_jones(summary, kernel)
  end
  throw(ArgumentError("unknown bandwidth rule :$rule (expected :normal or :sj)"))
end

end # module Bandwidth
