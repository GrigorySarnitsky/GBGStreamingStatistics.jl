# Bandwidth

Bandwidth selection. For the marginal (`Histogram1D`): the kernel- and degree-aware normal-reference rule
(the default; Silverman the Gaussian/degree-0 case), the data-driven Sheather–Jones (`rule = :sj` classic,
`:isj` improved-via-DCT), fit-based AIC (`rule = :aic`), and asymmetric-kernel Chen plug-in
(`rule = :plugin`) selectors (leave-one-out CV a recorded binned-ties gap). For the joint (`Histogram2D`):
the response-aware local plug-in ``h^\star(x)`` (`rule = :local`).

```@docs
bandwidth
```
