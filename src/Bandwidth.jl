# Bandwidth — §4.5: bandwidth selection (global rules + the response-aware local plug-in) and the
# post-hoc FFT smoothing helpers. The smoothing is a convolution of the binned summary with the
# kernel; the FFT computes it in O(B log B) and amortizes over bandwidth — transform the bin counts
# once, then each h is a multiply by the kernel's transform (`docs/ARCHITECTURE.md` §4.5).

module Bandwidth

using AbstractFFTs
using FFTW

export bandwidth


"""
    bandwidth(summary, channel; rule, x)

**Estimand.** the MSE-optimal bandwidth for the requested estimator/channel.

**Estimators.** global plug-in rules (Silverman; Sheather & Jones 1991) for a marginal density, and
the load-bearing **response-aware *local* plug-in** for conditional moments:

```math
h^\\star(x) = \\left[
  \\frac{R(K)\\,\\sigma^2(x)}{\\mu_2(K)^2\\, m''(x)^2\\, f_X(x)\\, n}
\\right]^{1/5},
```

which depends on the response's conditional variance ``\\sigma^2(x)`` **and** the regression
curvature ``m''(x)`` — so **each response channel gets its own ``h^\\star``**, all read back from the
*one* stored summary post hoc (this is why accumulation must be response-agnostic and rich). Here
``R(K) = \\int K^2`` is the kernel roughness and ``\\mu_2(K) = \\int u^2 K(u)\\,\\mathrm{d}u`` its
second moment.

The smoothing itself is post hoc and FFT-amortized over the bandwidth grid the local ``h^\\star(x)``
sweep needs; it is dwarfed by the per-sample accumulation hot loop, so it is a secondary performance
concern.

Reference: S. P. Sheather & M. C. Jones (1991), *A reliable data-based bandwidth selection method for
kernel density estimation*, J. R. Statist. Soc. B 53(3).

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function bandwidth end

end # module Bandwidth
