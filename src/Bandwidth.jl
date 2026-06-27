# Bandwidth — §4.5: bandwidth selection (global rules now; the response-aware local plug-in and the
# FFT-amortized smoothing are forthcoming). The smoothing is a convolution of the binned summary with
# the kernel; the FFT computes it in O(B log B) and amortizes over bandwidth (`docs/ARCHITECTURE.md`
# §4.5) — hence the FFTW dependency, used once that lands.

module Bandwidth

using AbstractFFTs
using FFTW

using ..Histograms: Histogram1D, _total
using ..DensityEstimation: gaussian_kernel, _normal_reference_h

export bandwidth


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

The data-driven **Sheather–Jones** plug-in (`rule = :sj`) and the fit-based **likelihood-CV / AIC**
selectors (`rule = :lcv`, `:aic`) are forthcoming (`CLAUDE.md` roadmap), as is the response-aware
*local* plug-in ``h^\\star(x) = \\bigl[R(K)\\,\\sigma^2(x) / (\\mu_2(K)^2\\, m''(x)^2\\, f_X(x)\\,
n)\\bigr]^{1/5}`` — each response channel its own bandwidth from the conditional variance
``\\sigma^2(x)`` and regression curvature ``m''(x)``.

References: B. W. Silverman (1986), [*Density Estimation for Statistics and Data
Analysis*](https://doi.org/10.1201/9781315140919), Chapman & Hall (the normal-reference rule and the
canonical bandwidth); C. R. Loader (1996), [*Local likelihood density
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
  end
  throw(ArgumentError("unknown bandwidth rule :$rule (expected :normal)"))
end

end # module Bandwidth
