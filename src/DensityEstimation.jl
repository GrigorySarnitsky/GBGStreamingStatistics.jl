# DensityEstimation — §4.1: the marginal density f_X, estimated from a 1-D histogram by binned KDE
# (WARPing). Written estimand-first (`docs/ARCHITECTURE.md` §2, §4.1).

module DensityEstimation

using ..Histograms: Histogram1D

export density


"""
    density(summary::Histogram1D, x; h)   # ≡ f_X

**Estimand.** the marginal (Lebesgue) density ``f_X = \\mathrm{d}\\mu_X/\\mathrm{d}\\lambda`` of the
conditioning variable ``X``, defined by ``P(X \\in A) = \\int_A f_X\\,\\mathrm{d}\\lambda``.

**Estimator** (binned KDE / WARPing). From the bin counts ``n_j`` at centers ``c_j`` of `summary`
(a [`Histogram1D`](@ref GBGStreamingStatistics.Histograms.Histogram1D)) and a kernel ``K`` with ``\\int K = 1`` and bandwidth ``h``,

```math
\\hat f_h(x) = \\frac{1}{n\\,h} \\sum_{j} n_j\\, K\\!\\left(\\frac{x - c_j}{h}\\right),
```

the binned approximation to ``\\frac{1}{nh}\\sum_i K\\!\\left(\\frac{x - X_i}{h}\\right)``. The
bandwidth ``h`` comes from [`bandwidth`](@ref GBGStreamingStatistics.Bandwidth.bandwidth); the smoothing is a convolution of the bin counts,
computed post hoc and amortized over ``h`` by the FFT.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function density end

end # module DensityEstimation
