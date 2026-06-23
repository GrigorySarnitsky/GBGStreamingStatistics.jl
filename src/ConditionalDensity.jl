# ConditionalDensity — §4.3: the conditional density f_{Y|X} and the two-time transition density
# p_τ(y|x), from a joint histogram by the double-kernel estimator (`docs/ARCHITECTURE.md` §4.3).

module ConditionalDensity

using ..Histograms: Histogram2D

export conditional_density, transition_density


"""
    conditional_density(summary::Histogram2D, x, y; h)   # ≡ f_{Y|X}

**Estimand.** the conditional density
``f_{Y\\mid X}(y \\mid x) = f_{X,Y}(x, y) / f_X(x)``.

**Estimator** (double-kernel; Rosenblatt 1969; Hyndman, Bashtannyk & Grunwald 1996). From the joint
histogram `summary` (a [`Histogram2D`](@ref GBGStreamingStatistics.Histograms.Histogram2D)), with bandwidths ``h_1`` (in ``x``) and ``h_2`` (in
``y``),

```math
\\hat f(y \\mid x) =
  \\frac{\\sum_i K_{h_1}(X_i - x)\\, K_{h_2}(Y_i - y)}{\\sum_i K_{h_1}(X_i - x)}.
```

The conditional MISE rate is slower than a marginal density (``O(n^{-2/3})``).

Reference: M. Rosenblatt (1969), *Conditional probability density and regression estimators*, in
*Multivariate Analysis II*.

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function conditional_density end


"""
    transition_density(summary::Histogram2D, x, y)   # ≡ p_τ(y|x)

**Estimand.** the **two-time transition density**
``p_\\tau(y \\mid x) = f_{X_{t+\\tau} \\mid X_t}(y \\mid x)`` — the conditional density of the state a
lag ``\\tau`` later, given the present state (nonstationary: ``p_{t,\\tau}``). For the
Ornstein–Uhlenbeck oracle this is known in closed form (the natural end-to-end test).

**Estimator.** the [`conditional_density`](@ref GBGStreamingStatistics.ConditionalDensity.conditional_density) double-kernel estimator applied to the per-lag joint
histogram of pairs ``\\bigl(X_t,\\, X_{t+\\tau}\\bigr)`` (binned by [`observe_pairs!`](@ref GBGStreamingStatistics.Histograms.observe_pairs!)).

!!! note "Scaffold stage"
    Declared; not yet implemented.
"""
function transition_density end

end # module ConditionalDensity
