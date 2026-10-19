# Bandwidth

Bandwidths are selected from a summary, the binned histogram of the samples. A marginal
[`Histogram1D`](@ref) supports, for [`density`](@ref), the kernel- and degree-aware normal-reference rule
(the default; Silverman's is its Gaussian, degree-0 case), Sheather–Jones (`rule = :sj`; `:isj`
improved via the discrete cosine transform; preferred), fit-based AIC (`rule = :aic`), and the
asymmetric-kernel Chen plug-in (`rule = :plugin`). A joint [`Histogram2D`](@ref) supports the
response-aware local plug-in ``h^\star(x)`` (`rule = :local`), one bandwidth per conditioning value.

```@docs
bandwidth
```
