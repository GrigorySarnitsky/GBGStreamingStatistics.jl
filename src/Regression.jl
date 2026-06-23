# Regression — §4.2: the conditional mean m(x) and conditional variance σ²(x), from per-bin
# power-sums by local-linear regression. The estimand-first reference module — write every estimator
# docstring to this pattern (`docs/ARCHITECTURE.md` §2, §4.2).

module Regression

using LinearAlgebra

using ..Histograms: Histogram1D

export conditional_mean, conditional_variance


"""
    conditional_mean(summary::Histogram1D, x; h)   # ≡ m(x)

**Estimand.** the conditional expectation
``m(x) = \\mathbb{E}[Y \\mid X = x] = \\int y\\, f_{Y\\mid X}(y \\mid x)\\,\\mathrm{d}y`` — the
``L^2(\\sigma(X))`` projection of the response ``Y`` onto the conditioning variable ``X``, i.e. the
regression function.

**Estimator** (local-linear; Fan & Gijbels). From the per-bin power-sums of `summary` form the
kernel-weighted sums ``T_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell`` and
``P_\\ell(x) = \\sum_i K_h(X_i-x)(X_i-x)^\\ell Y_i`` and solve the local normal equations:

```math
\\hat m_h(x) = \\frac{T_2 P_0 - T_1 P_1}{T_0 T_2 - T_1^{2}}.
```

Nadaraya–Watson is the ``T_1 = T_2 = 0`` reduction ``\\hat m_h = P_0/T_0`` (kept for reference;
local-linear is the default — it removes the design-density bias ``2 m'(x) f'(x)/f(x)``). The
consistency claim (`docs/FORMALIZATION.md`):
``\\operatorname{bias}\\hat m_h(x) \\approx \\tfrac12 h^2 \\mu_2(K)\\, m''(x)`` and
``\\operatorname{Var}\\hat m_h(x) \\approx R(K)\\,\\sigma^2(x)/(n\\,h\\,f_X(x))``.

Reference: J. Fan & I. Gijbels (1996), *Local Polynomial Modelling and Its Applications*,
Chapman & Hall.

!!! note "Forthcoming"
    Declared; not yet implemented.
"""
function conditional_mean end


"""
    conditional_variance(summary::Histogram1D, x; h)   # ≡ σ²(x)

**Estimand.** the conditional variance
``\\sigma^2(x) = \\operatorname{Var}(Y \\mid X = x) = \\mathbb{E}[Y^2 \\mid X = x] - m(x)^2``, with
``m(x)`` the [`conditional_mean`](@ref GBGStreamingStatistics.Regression.conditional_mean).

**Estimator.** the same local-linear fit applied to the response ``Y^2``, combined as
``\\widehat{\\mathbb{E}}[Y^2\\mid x] - \\hat m_h(x)^2`` (Fan & Yao 1998). For heavy-tailed channels
(e.g. ``(\\Delta X)^6``) a one-step robust local M-estimator (Fan & Jiang 2000) replaces the
least-squares fit.

Reference: J. Fan & Q. Yao (1998), *Efficient estimation of conditional variance functions in
stochastic regression*, Biometrika 85(3).

!!! note "Forthcoming"
    Declared; not yet implemented.
"""
function conditional_variance end

end # module Regression
