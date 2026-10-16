# Bandwidth

Bandwidth selection for the [`density`](@ref) estimators. For the marginal
([`Histogram1D`](Histograms.md)): the kernel- and degree-aware normal-reference rule (the default;
Silverman's rule is its Gaussian, degree-0 case), the data-driven Sheather–Jones (`rule = :sj` classic,
`:isj` improved via the discrete cosine transform), fit-based AIC (`rule = :aic`), and asymmetric-kernel Chen plug-in (`rule = :plugin`)
selectors; for a data-driven choice prefer `rule = :sj` or `:isj`. For the joint (`Histogram2D`): the
response-aware local plug-in ``h^\star(x)`` (`rule = :local`).

```@docs
bandwidth
```
