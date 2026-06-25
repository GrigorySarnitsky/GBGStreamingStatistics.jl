# Bandwidth — §4.5: bandwidth selection (global rules now; the response-aware local plug-in and the
# FFT-amortized smoothing are forthcoming). The smoothing is a convolution of the binned summary with
# the kernel; the FFT computes it in O(B log B) and amortizes over bandwidth (`docs/ARCHITECTURE.md`
# §4.5) — hence the FFTW dependency, used once that lands.

module Bandwidth

using AbstractFFTs
using FFTW

using ..Histograms: Histogram1D, _total, _sample_std

export bandwidth


"""
    bandwidth(summary::Histogram1D; rule = :silverman) -> h

**Estimand.** the MSE-optimal smoothing bandwidth for the marginal density of `summary`.

**Estimator** (global plug-in rule of thumb). From the recovered standard deviation ``\\hat\\sigma``
of the binned variable and the sample size ``n``,

```math
h = c_{\\text{rule}}\\,\\hat\\sigma\\,n^{-1/5},
\\qquad
c_{\\texttt{silverman}} = (4/3)^{1/5} \\approx 1.0592,\\quad c_{\\texttt{scott}} = 1,
```

the AMISE-optimal bandwidth for a Gaussian kernel under a Gaussian reference (`rule = :silverman`,
the default; `:scott` drops the constant). Both read ``\\hat\\sigma`` and ``n`` back from the *one*
stored summary, post hoc.

The **response-aware local plug-in** ``h^\\star(x) = \\bigl[R(K)\\,\\sigma^2(x) / (\\mu_2(K)^2\\,
m''(x)^2\\, f_X(x)\\, n)\\bigr]^{1/5}`` — which gives each response channel its own bandwidth from the
conditional variance ``\\sigma^2(x)`` and regression curvature ``m''(x)`` — is forthcoming with the
conditional estimators (`CLAUDE.md` roadmap); here ``R(K) = \\int K^2`` and
``\\mu_2(K) = \\int u^2 K(u)\\,\\mathrm{d}u``.

Reference: B. W. Silverman (1986), [*Density Estimation for Statistics and Data Analysis*](https://doi.org/10.1201/9781315140919), Chapman &
Hall.
"""
function bandwidth(summary::Histogram1D; rule::Symbol = :silverman)
  n = _total(summary)
  n ≥ 2 || throw(ArgumentError("bandwidth needs at least 2 samples (got n = $(Int(n)))"))
  c = rule === :silverman ? (4 / 3)^(1 / 5) :
      rule === :scott     ? 1.0 :
      throw(ArgumentError("unknown bandwidth rule :$rule (expected :silverman or :scott)"))
  return c * _sample_std(summary) * float(n)^(-1 / 5)
end

end # module Bandwidth
